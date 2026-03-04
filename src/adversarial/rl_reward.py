"""RL reward shaping for the Layer 7 Bayesian Graph VAE.

Implements a REINFORCE (policy-gradient) reward module that shapes the
VAE training objective using:

    reward(molecule) = binding_affinity_score(L8) − constraint_penalty(L9)

where:
  • binding_affinity_score comes from the Layer 8 PropertyPredictor's
    QED prediction (the best drug-likeness proxy available from L8).
  • constraint_penalty comes from the ConstraintViolationVector built
    from the Rocq bridge result (Step 2).

Gradient path
─────────────
The Rocq proof checker is a black-box oracle — no gradients flow through
it.  REINFORCE uses the log-probability of the VAE's latent sample z
under the encoder distribution q(z|x) as the surrogate loss:

    L_RL = −(reward − baseline) · log q(z | μ, log σ²)

where:
    log q(z | μ, log σ²) = −½ · Σ_d [(z_d − μ_d)² / σ_d² + log σ_d² + log 2π]

The baseline is an exponential moving average of recent rewards, which
reduces gradient variance without introducing bias (Williams 1992;
Olivecrona et al. 2017 / REINVENT).

Integration with model.py
─────────────────────────
The existing vae_loss_function is augmented:

    L_total = L_VAE(recon, KL) + λ_RL · L_RL

See model.py → train_with_rl_reward() for the complete training loop.
"""

import math
import os
from typing import Optional

import torch
import torch.nn as nn

from .rocq_bridge import check_molecule
from .constraint_vector import ConstraintViolationVector, load_weights


# ---------------------------------------------------------------------------
# Log-probability of z under the encoder distribution
# ---------------------------------------------------------------------------

def log_q_gaussian(z: torch.Tensor, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
    """Per-sample log-probability log q(z | μ, log σ²).

    Parameters
    ----------
    z, mu, logvar : (batch, latent_dim)

    Returns
    -------
    log_q : (batch,)  — sum over latent dimensions
    """
    return -0.5 * (
        (z - mu).pow(2) / logvar.exp()
        + logvar
        + math.log(2 * math.pi)
    ).sum(dim=1)


# ---------------------------------------------------------------------------
# REINFORCE loss
# ---------------------------------------------------------------------------

def reinforce_loss(
    log_q: torch.Tensor,
    rewards: torch.Tensor,
    baseline: float = 0.0,
) -> torch.Tensor:
    """REINFORCE policy-gradient loss (Williams 1992).

    L = −E[(reward − baseline) · log q(z | x)]

    Parameters
    ----------
    log_q   : (batch,)  encoder log-probabilities from log_q_gaussian()
    rewards : (batch,)  per-molecule reward values
    baseline: float     running mean reward subtracted for variance reduction

    Returns
    -------
    Scalar loss tensor (mean over batch).
    """
    advantages = rewards.detach() - baseline  # shape: (batch,)
    return -(advantages * log_q).mean()


# ---------------------------------------------------------------------------
# RewardShaper
# ---------------------------------------------------------------------------

class RewardShaper:
    """Computes per-molecule rewards and the REINFORCE loss term.

    Parameters
    ----------
    weights_config : str, optional
        Path to constraint_weights.yaml.  Defaults to the repo-level config.
    baseline_ema_alpha : float
        Smoothing factor for the exponential moving average baseline.
    property_predictor : nn.Module, optional
        Trained Layer 8 PropertyPredictor for QED scoring.  If None, the
        QED contribution is skipped (pure constraint-penalty reward).
    """

    def __init__(
        self,
        weights_config: Optional[str] = None,
        baseline_ema_alpha: float = 0.05,
        property_predictor: Optional[nn.Module] = None,
    ):
        self.weights = load_weights(weights_config)
        self.baseline = 0.0
        self._alpha = baseline_ema_alpha
        self._property_predictor = property_predictor

        # Scaler tensors for denormalising PropertyPredictor output
        self._scaler_mean: Optional[torch.Tensor] = None
        self._scaler_std:  Optional[torch.Tensor] = None

    # ------------------------------------------------------------------
    # Reward computation
    # ------------------------------------------------------------------

    def reward_for_smiles(self, smiles: str, binding_affinity: float = 0.0) -> float:
        """Compute the scalar reward for a single molecule.

        Parameters
        ----------
        smiles           : SMILES string of the molecule.
        binding_affinity : Pre-computed binding affinity score (e.g. QED
                           from PropertyPredictor).  Pass 0.0 to use only
                           the constraint penalty.

        Returns
        -------
        reward = binding_affinity − constraint_penalty
        """
        result = check_molecule(smiles)
        vec = ConstraintViolationVector.from_bridge_result(result)
        penalty = vec.penalty_score(self.weights)
        return binding_affinity - penalty

    def compute_batch_rewards(
        self,
        smiles_batch: list,
        fingerprints: Optional[torch.Tensor] = None,
    ) -> tuple:
        """Compute rewards for a batch of SMILES strings.

        If a PropertyPredictor is available and fingerprints are supplied,
        the QED prediction is used as the binding_affinity term.

        Parameters
        ----------
        smiles_batch  : list[str]  — one SMILES per molecule
        fingerprints  : (batch, 2048) tensor, optional

        Returns
        -------
        rewards           : (batch,) float tensor
        violation_rates   : dict[str, float]  — fraction violated per type
        """
        # --- Optional QED from PropertyPredictor ---
        qed_scores = self._predict_qed(fingerprints, len(smiles_batch))

        rewards_list = []
        violation_counts = {c: 0 for c in self.weights}

        for i, smiles in enumerate(smiles_batch):
            result = check_molecule(smiles)
            vec = ConstraintViolationVector.from_bridge_result(result)
            penalty = vec.penalty_score(self.weights)
            r = float(qed_scores[i]) - penalty
            rewards_list.append(r)

            for ctype, viol in vec.violated.items():
                if viol and ctype in violation_counts:
                    violation_counts[ctype] += 1

        n = max(len(smiles_batch), 1)
        violation_rates = {c: violation_counts[c] / n for c in violation_counts}

        return torch.tensor(rewards_list, dtype=torch.float32), violation_rates

    def update_baseline(self, mean_reward: float) -> None:
        """Update the EMA baseline with the mean reward of the last batch."""
        self.baseline = (1 - self._alpha) * self.baseline + self._alpha * mean_reward

    def load_property_predictor(self, ckpt_path: str) -> None:
        """Load a saved PropertyPredictor checkpoint (from src/ml/property_predictor.py)."""
        import sys
        repo_root = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..")
        )
        if repo_root not in sys.path:
            sys.path.insert(0, repo_root)

        from src.ml.property_predictor import PropertyPredictor

        ckpt = torch.load(ckpt_path, weights_only=True)
        model = PropertyPredictor()
        model.load_state_dict(ckpt["model"])
        model.eval()

        self._property_predictor = model
        self._scaler_mean = ckpt.get("scaler_mean")
        self._scaler_std  = ckpt.get("scaler_std")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _predict_qed(
        self,
        fingerprints: Optional[torch.Tensor],
        n: int,
    ) -> list:
        """Return a list of QED scores, one per molecule.

        QED is index 0 of PROPERTY_NAMES = ["QED", "LogP", "MW", "TPSA", "HBD", "HBA"].
        Falls back to 0.5 (neutral) when the predictor is unavailable.
        """
        if self._property_predictor is None or fingerprints is None:
            return [0.5] * n

        self._property_predictor.train()  # keep dropout for MC sampling
        with torch.no_grad():
            preds = self._property_predictor(fingerprints)  # (batch, 6)

        # Denormalise if scaler tensors are available
        if self._scaler_mean is not None and self._scaler_std is not None:
            preds = preds * (self._scaler_std + 1e-8) + self._scaler_mean

        qed = preds[:, 0].tolist()  # QED is the first property
        # Clamp to [0, 1] since QED is a drug-likeness score in that range
        return [max(0.0, min(1.0, q)) for q in qed]
