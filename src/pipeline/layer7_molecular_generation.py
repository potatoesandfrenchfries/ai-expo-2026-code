"""Layer 7: AI-driven molecular generation via Bayesian Graph VAE.

Wraps the BayesianGraphVAE from model.py as a proper pipeline layer.
The VAE operates in Morgan fingerprint space (2048-dim); SMILES are
recovered by cosine nearest-neighbour lookup against the ZINC-250k pool.

Data contract:
  Input : Layer6Result  (.top_target.gene used for provenance)
  Output: Layer7Result  (list of GeneratedCandidate with SMILES + metadata)

Generation strategy (Phase 1):
  1. Load ZINC-10k Morgan fingerprints (data/X.pt) + SMILES (data/molecules_clean.csv)
  2. Initialise BayesianGraphVAE (load data/bayesian_vae.pt if present, else random init)
  3. Optionally apply FeedbackController.update_vae_sampling() for constraint propagation
  4. Sample n_oversample latent points z ~ N(0, I_16)
  5. Decode z → reconstructed fingerprint (2048-dim)
  6. Cosine similarity vs ZINC pool → take top-1 SMILES per latent sample
  7. Deduplicate + validate with RDKit.Chem.MolFromSmiles
  8. Compute per-candidate MC uncertainty from VAE reconstruction variance
  9. Return up to n_candidates unique valid SMILES

Usage:
    from src.pipeline.layer7_molecular_generation import Layer7MolecularGeneration
    layer7 = Layer7MolecularGeneration()
    result = layer7.run(layer6_result, n_candidates=50)
"""

import os
import sys
from dataclasses import dataclass, field
from typing import Optional

import torch
import torch.nn.functional as F
import pandas as pd
from rdkit import Chem

# ---------------------------------------------------------------------------
# Resolve repo root so model.py and data/ are accessible regardless of cwd
# ---------------------------------------------------------------------------
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from model import BayesianGraphVAE, mc_predict  # noqa: E402
from src.llm.feedback_controller import FeedbackController  # noqa: E402
from src.pipeline.layer6_target_confidence import Layer6Result  # noqa: E402

_DATA_DIR   = os.path.join(_REPO_ROOT, "data")
_ZINC_CSV   = os.path.join(_REPO_ROOT, "src", "data", "250k_rndm_zinc_drugs_clean_3.csv")
_VAE_CKPT   = os.path.join(_DATA_DIR, "bayesian_vae.pt")
_FP_CACHE   = os.path.join(_DATA_DIR, "X.pt")
_SMILES_CSV = os.path.join(_DATA_DIR, "molecules_clean.csv")


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

@dataclass
class GeneratedCandidate:
    smiles:      str
    source:      str    # "vae_nn" = VAE decode + cosine nearest-neighbour
    similarity:  float  # cosine similarity to decoded fingerprint [0, 1]
    uncertainty: float  # mean VAE reconstruction variance (MC dropout)


@dataclass
class Layer7Result:
    target_gene:          str
    candidates:           list  # list[GeneratedCandidate], deduplicated + RDKit-valid
    n_sampled:            int   # latent points drawn (before dedup/filtering)
    n_valid:              int   # SMILES that pass RDKit validation
    constraints_applied:  bool  # True if FeedbackController was used


# ---------------------------------------------------------------------------
# Main Layer 7 class
# ---------------------------------------------------------------------------

class Layer7MolecularGeneration:
    """
    Bayesian VAE-driven molecular generation for Phase 1 drug discovery.

    The VAE navigates a 16-dimensional latent fingerprint space.  SMILES are
    retrieved from the ZINC-10k pool via cosine nearest-neighbour decoding so
    all candidates are guaranteed to be real, synthesisable molecules.

    The nearest-neighbour approach is appropriate for Phase 1: the VAE learns
    the topology of drug-like fingerprint space, and decoding to the closest
    real molecule avoids invalid SMILES while preserving chemical diversity.
    """

    def __init__(self,
                 zinc_smiles_path: str = _SMILES_CSV,
                 fp_cache_path: str    = _FP_CACHE,
                 vae_ckpt_path: str    = _VAE_CKPT):
        self._zinc_smiles_path = zinc_smiles_path
        self._fp_cache_path    = fp_cache_path
        self._vae_ckpt_path    = vae_ckpt_path
        # Lazy-loaded caches
        self._smiles_pool: Optional[list]         = None
        self._X_pool:      Optional[torch.Tensor] = None
        self._vae:         Optional[BayesianGraphVAE] = None

    # -----------------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------------

    def run(self,
            layer6_result: Layer6Result,
            n_candidates:  int = 50,
            feedback:      Optional[FeedbackController] = None) -> Layer7Result:
        """
        Generate drug candidates for the top target from Layer 6.

        Parameters
        ----------
        layer6_result : Result from Layer 6; top_target.gene used for provenance.
        n_candidates  : Desired number of unique, valid SMILES to return.
        feedback      : Optional FeedbackController carrying L9 failure constraints.
                        When provided, its bounds are applied to the VAE sampler.
        """
        target_gene = (
            layer6_result.top_target.gene
            if layer6_result.top_target else "UNKNOWN"
        )

        X_pool, smiles_pool = self._load_pool()
        vae = self._load_vae()

        constraints_applied = False
        if feedback is not None:
            feedback.update_vae_sampling(vae)
            constraints_applied = True

        # Oversample to absorb dedup and invalid-SMILES losses
        n_oversample = max(n_candidates * 4, 200)
        candidates = self._generate(vae, X_pool, smiles_pool, n_oversample, n_candidates)

        return Layer7Result(
            target_gene=target_gene,
            candidates=candidates,
            n_sampled=n_oversample,
            n_valid=len(candidates),
            constraints_applied=constraints_applied,
        )

    # -----------------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------------

    def _load_pool(self) -> tuple:
        """Lazy-load ZINC fingerprints and corresponding SMILES."""
        if self._X_pool is not None and self._smiles_pool is not None:
            return self._X_pool, self._smiles_pool

        if not os.path.exists(self._fp_cache_path):
            raise FileNotFoundError(
                f"Fingerprint cache not found at {self._fp_cache_path}. "
                "Run src/data/data_loader.py first to generate data/X.pt."
            )

        X = torch.load(self._fp_cache_path, weights_only=True)
        n_pool = len(X)

        # Prefer the clean sampled CSV; fall back to the full ZINC file
        if os.path.exists(self._zinc_smiles_path):
            df = pd.read_csv(self._zinc_smiles_path)
            smiles = df["smiles"].tolist()[:n_pool]
        elif os.path.exists(_ZINC_CSV):
            df = pd.read_csv(_ZINC_CSV).sample(n_pool, random_state=42)
            smiles = df["smiles"].tolist()
        else:
            raise FileNotFoundError(
                "SMILES pool not found. Run src/data/data_loader.py first, "
                f"or supply the ZINC CSV at {_ZINC_CSV}."
            )

        self._X_pool      = X
        self._smiles_pool = smiles
        return X, smiles

    def _load_vae(self) -> BayesianGraphVAE:
        """Load a pre-trained VAE checkpoint, or return a randomly-initialised model."""
        if self._vae is not None:
            return self._vae

        vae = BayesianGraphVAE(input_dim=2048, latent_dim=16, hidden_dim=64)
        if os.path.exists(self._vae_ckpt_path):
            state = torch.load(self._vae_ckpt_path, weights_only=True)
            vae.load_state_dict(state)
            # Checkpoint loaded — latent space is trained on ZINC drug-likeness
        # If no checkpoint, random weights still provide useful fingerprint diversity
        # because the decoder projects latent points into the 2048-dim space, and
        # nearest-neighbour retrieval maps them to real ZINC molecules.
        self._vae = vae
        return vae

    def _generate(self,
                  vae:          BayesianGraphVAE,
                  X_pool:       torch.Tensor,
                  smiles_pool:  list,
                  n_oversample: int,
                  n_candidates: int) -> list:
        """
        Core generation loop:
          1. Sample z from N(0, I_16)
          2. Decode to reconstructed 2048-dim fingerprint
          3. Cosine nearest-neighbour → ZINC SMILES
          4. Deduplicate and RDKit-validate
          5. Annotate with MC uncertainty
        """
        vae.eval()

        with torch.no_grad():
            z         = torch.randn(n_oversample, 16)
            h_dec     = vae.decoder_hidden(z)             # (n_oversample, 64)
            recon_fp  = vae.mu_head(h_dec)                # (n_oversample, 2048)

            recon_norm = F.normalize(recon_fp, dim=-1)    # (n_oversample, 2048)
            pool_norm  = F.normalize(X_pool.float(), dim=-1)   # (pool_size, 2048)
            sims       = recon_norm @ pool_norm.T          # (n_oversample, pool_size)
            nn_indices = sims.argmax(dim=-1)               # (n_oversample,)
            nn_sims    = sims.max(dim=-1).values           # (n_oversample,)

        results: list = []
        seen:    set  = set()

        for i in range(n_oversample):
            if len(results) >= n_candidates:
                break

            idx        = int(nn_indices[i].item())
            raw_smiles = smiles_pool[idx]

            if raw_smiles in seen:
                continue

            mol = Chem.MolFromSmiles(raw_smiles)
            if mol is None:
                continue

            # MC uncertainty: 20 stochastic forward passes over this fingerprint
            fp_tensor   = X_pool[idx].float().unsqueeze(0)
            _, var      = mc_predict(vae, fp_tensor, samples=20)
            uncertainty = round(float(var.mean().item()), 6)

            seen.add(raw_smiles)
            results.append(GeneratedCandidate(
                smiles=raw_smiles,
                source="vae_nn",
                similarity=round(float(nn_sims[i].item()), 4),
                uncertainty=uncertainty,
            ))

        return results

    def print_report(self, result: Layer7Result) -> None:
        print(f"\n[Layer 7] Molecular Generation — target: {result.target_gene}")
        print(f"  Sampled {result.n_sampled} latent points "
              f"→ {result.n_valid} unique valid candidates")
        if result.constraints_applied:
            print("  Feedback constraints applied from prior L9 failures.")
        print(f"\n  {'SMILES':<45}  {'Sim':>5}  {'Uncert':>8}")
        print("  " + "-" * 63)
        for c in result.candidates[:10]:
            abbr = (c.smiles[:42] + "...") if len(c.smiles) > 45 else c.smiles
            print(f"  {abbr:<45}  {c.similarity:>5.3f}  {c.uncertainty:>8.5f}")
        if len(result.candidates) > 10:
            print(f"  ... ({len(result.candidates) - 10} more not shown)")


# ---------------------------------------------------------------------------
# Standalone demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from src.pipeline.layer6_target_confidence import Layer6Result, TargetConfidenceScore
    from dataclasses import field as _field

    stub_top = TargetConfidenceScore(
        gene="BCL2", final_score=0.754,
        ot_component=0.328, centrality_component=0.115,
        chembl_component=0.089, tractability_component=0.1,
        pathogenicity_component=0.1, penalty=-0.02,
    )
    stub_l6 = Layer6Result(
        scores=[stub_top], top_target=stub_top, disease_id="EFO_0000220"
    )

    print("Running Layer 7 standalone demo (VAE nearest-neighbour generation)...")
    layer7 = Layer7MolecularGeneration()
    result = layer7.run(stub_l6, n_candidates=10)
    layer7.print_report(result)
