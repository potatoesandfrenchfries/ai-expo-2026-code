"""Surrogate discriminator — fast Rocq pass/fail pre-filter (Step 4).

Because a single Rocq proof check takes ~seconds and the generator may
produce thousands of candidates per training step, a lightweight MLP
surrogate can be trained to predict the Rocq verdict from a 2048-bit
Morgan fingerprint.

Workflow
────────
1. Accumulate (fingerprint, rocq_pass) labels from the Rocq bridge.
2. Periodically retrain the surrogate on all accumulated labels.
3. Use the surrogate to pre-score new candidates:
      fast_score = surrogate.predict_proba(fp)
4. Only send the top-k% (highest surrogate pass-probability) candidates
   to the full Rocq checker, saving expensive subprocess time.
5. Track surrogate accuracy vs. ground-truth; alert if it drops below
   the configured threshold.

Usage
─────
    disc = SurrogateDiscriminator()
    disc.add_labels(fingerprints, labels)     # accumulate Rocq verdicts
    disc.retrain()                            # fit on all accumulated data

    # Pre-filter: keep only top-50% by predicted pass probability
    fp_tensor = ...                           # (N, 2048)
    mask = disc.top_k_mask(fp_tensor, k=0.5) # bool tensor (N,)

    # Evaluate accuracy against held-out Rocq verdicts
    acc = disc.accuracy(fp_tensor, true_labels)
"""

import os
from collections import deque
from typing import Optional

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim


# ---------------------------------------------------------------------------
# MLP architecture
# ---------------------------------------------------------------------------

class _SurrogateMLP(nn.Module):
    """2-layer MLP binary classifier on 2048-dim Morgan fingerprints."""

    def __init__(self, input_dim: int = 2048, hidden_dim: int = 256, dropout: float = 0.2):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Returns raw logits (batch,)."""
        return self.net(x).squeeze(1)


# ---------------------------------------------------------------------------
# SurrogateDiscriminator
# ---------------------------------------------------------------------------

class SurrogateDiscriminator:
    """Train-and-predict surrogate for Rocq pass/fail labels.

    Parameters
    ----------
    input_dim          : fingerprint dimension (default 2048, Morgan r=2)
    hidden_dim         : MLP hidden layer width
    max_buffer_size    : maximum number of labelled examples to keep
    accuracy_threshold : alert if surrogate accuracy drops below this
    retrain_every_n    : retrain after this many new labels are added
    """

    def __init__(
        self,
        input_dim: int = 2048,
        hidden_dim: int = 256,
        max_buffer_size: int = 50_000,
        accuracy_threshold: float = 0.85,
        retrain_every_n: int = 500,
    ):
        self._model = _SurrogateMLP(input_dim, hidden_dim)
        self._trained = False

        # Replay buffer: deque of (fp_np_array, label) pairs
        self._buffer: deque = deque(maxlen=max_buffer_size)
        self._labels_since_last_train = 0

        self.accuracy_threshold = accuracy_threshold
        self.retrain_every_n = retrain_every_n

        # Running accuracy estimate (updated after each retrain)
        self.last_accuracy: float = float("nan")

    # ------------------------------------------------------------------
    # Label accumulation
    # ------------------------------------------------------------------

    def add_labels(
        self,
        fingerprints: torch.Tensor,
        labels: torch.Tensor,
    ) -> None:
        """Add (fingerprint, rocq_pass) pairs to the replay buffer.

        Parameters
        ----------
        fingerprints : (N, 2048) float tensor — Morgan fingerprints
        labels       : (N,)      float/bool tensor — 1.0 = passed, 0.0 = failed
        """
        fp_np = fingerprints.detach().cpu().numpy()
        lab_np = labels.detach().cpu().float().numpy()
        for i in range(len(fp_np)):
            self._buffer.append((fp_np[i], lab_np[i]))
        self._labels_since_last_train += len(fp_np)

    def add_smiles_labels(
        self,
        smiles_list: list,
        labels: list,
    ) -> None:
        """Convenience method: compute fingerprints from SMILES then add.

        Parameters
        ----------
        smiles_list : list[str]
        labels      : list[bool | int]  — True/1 = Rocq passed
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem, DataStructs
        fps = []
        valid_labels = []
        for smi, lab in zip(smiles_list, labels):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            arr = np.zeros(2048, dtype=np.float32)
            DataStructs.ConvertToNumpyArray(fp, arr)
            fps.append(arr)
            valid_labels.append(float(lab))

        if fps:
            fp_tensor = torch.tensor(np.stack(fps))
            lab_tensor = torch.tensor(valid_labels)
            self.add_labels(fp_tensor, lab_tensor)

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    def should_retrain(self, step: Optional[int] = None) -> bool:
        """Return True if enough new labels have accumulated to warrant a retrain."""
        if len(self._buffer) < 32:
            return False
        if step is not None:
            return step % self.retrain_every_n == 0
        return self._labels_since_last_train >= self.retrain_every_n

    def retrain(
        self,
        epochs: int = 20,
        batch_size: int = 256,
        lr: float = 1e-3,
        validation_fraction: float = 0.1,
    ) -> float:
        """Fit the surrogate on the current replay buffer.

        Returns validation accuracy (fraction correctly classified).
        """
        if len(self._buffer) < 32:
            return float("nan")

        fps = np.stack([b[0] for b in self._buffer])
        labs = np.array([b[1] for b in self._buffer], dtype=np.float32)

        n = len(fps)
        n_val = max(16, int(n * validation_fraction))
        idx = np.random.permutation(n)
        val_idx, train_idx = idx[:n_val], idx[n_val:]

        X_train = torch.tensor(fps[train_idx])
        y_train = torch.tensor(labs[train_idx])
        X_val   = torch.tensor(fps[val_idx])
        y_val   = torch.tensor(labs[val_idx])

        self._model.train()
        optimizer = optim.Adam(self._model.parameters(), lr=lr)
        criterion = nn.BCEWithLogitsLoss()

        for _ in range(epochs):
            perm = torch.randperm(len(X_train))
            for i in range(0, len(X_train), batch_size):
                xb = X_train[perm[i: i + batch_size]]
                yb = y_train[perm[i: i + batch_size]]
                optimizer.zero_grad()
                logits = self._model(xb)
                loss = criterion(logits, yb)
                loss.backward()
                optimizer.step()

        # Validation accuracy
        self._model.eval()
        with torch.no_grad():
            logits_val = self._model(X_val)
            preds = (torch.sigmoid(logits_val) > 0.5).float()
            acc = (preds == y_val).float().mean().item()

        self.last_accuracy = acc
        self._trained = True
        self._labels_since_last_train = 0

        if acc < self.accuracy_threshold:
            print(
                f"[SurrogateDiscriminator] WARNING: validation accuracy "
                f"{acc:.3f} < threshold {self.accuracy_threshold:.3f}. "
                f"Consider increasing the training corpus or retraining more frequently."
            )

        return acc

    # ------------------------------------------------------------------
    # Inference
    # ------------------------------------------------------------------

    def predict_proba(self, fingerprints: torch.Tensor) -> torch.Tensor:
        """Return Rocq-pass probabilities for a batch of fingerprints.

        Returns
        -------
        (batch,) float tensor in [0, 1].  Returns 0.5 (neutral) if not
        yet trained.
        """
        if not self._trained:
            return torch.full((len(fingerprints),), 0.5)

        self._model.eval()
        with torch.no_grad():
            logits = self._model(fingerprints)
            return torch.sigmoid(logits)

    def top_k_mask(self, fingerprints: torch.Tensor, k: float = 0.5) -> torch.Tensor:
        """Boolean mask of the top-k% candidates by predicted pass probability.

        Parameters
        ----------
        fingerprints : (N, 2048)
        k            : fraction to keep (e.g. 0.5 keeps the top 50%)

        Returns
        -------
        (N,) bool tensor — True for candidates that should proceed to
        full Rocq checking.
        """
        proba = self.predict_proba(fingerprints)
        if not self._trained:
            return torch.ones(len(fingerprints), dtype=torch.bool)
        n_keep = max(1, int(len(proba) * k))
        threshold = proba.topk(n_keep).values.min().item()
        return proba >= threshold

    def accuracy(
        self,
        fingerprints: torch.Tensor,
        true_labels: torch.Tensor,
        threshold: float = 0.5,
    ) -> float:
        """Compute binary classification accuracy against ground-truth labels.

        Parameters
        ----------
        fingerprints : (N, 2048)
        true_labels  : (N,) bool/float — 1.0 = Rocq passed
        threshold    : sigmoid threshold for positive class

        Returns
        -------
        Fraction of correctly classified molecules.
        """
        proba = self.predict_proba(fingerprints)
        preds = (proba > threshold).float()
        true  = true_labels.float()
        return (preds == true).float().mean().item()

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: str) -> None:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        torch.save(
            {
                "model_state": self._model.state_dict(),
                "trained": self._trained,
                "last_accuracy": self.last_accuracy,
            },
            path,
        )

    def load(self, path: str) -> None:
        ckpt = torch.load(path, weights_only=True)
        self._model.load_state_dict(ckpt["model_state"])
        self._trained = ckpt.get("trained", True)
        self.last_accuracy = ckpt.get("last_accuracy", float("nan"))
