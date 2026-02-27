"""Layer 8: Molecular property prediction.

Multi-output MLP trained on 2048-dim Morgan fingerprints (from data_loader.py)
to predict six physicochemical properties simultaneously:
  • QED        — drug-likeness score [0, 1]
  • LogP       — lipophilicity
  • MW         — molecular weight (Da)
  • TPSA       — topological polar surface area (Å²)
  • HBD        — H-bond donor count
  • HBA        — H-bond acceptor count

All properties are predicted with an independent regression head so they can
be optimised or constrained separately.  Epistemic uncertainty is estimated
via MC Dropout (model.train() at inference time).

Usage:
    python -m src.ml.property_predictor
"""

import os
from dataclasses import dataclass

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA_DIR  = os.path.join(REPO_ROOT, "data")

PROPERTY_NAMES = ["QED", "LogP", "MW", "TPSA", "HBD", "HBA"]
N_PROPS = len(PROPERTY_NAMES)


# ---------------------------------------------------------------------------
# Label builder  (runs once; reads molecules_clean.csv → computes targets)
# ---------------------------------------------------------------------------

def build_property_labels(csv_path: str) -> torch.Tensor:
    """
    Compute a (N, 6) tensor of property labels from molecules_clean.csv.
    Column order matches PROPERTY_NAMES.
    """
    import pandas as pd

    df = pd.read_csv(csv_path)
    rows = []
    for smi in df["smiles"]:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rows.append([0.0] * N_PROPS)
            continue
        rows.append([
            float(df.loc[df["smiles"] == smi, "qed"].iloc[0]),
            Descriptors.MolLogP(mol),
            Descriptors.ExactMolWt(mol),
            Descriptors.TPSA(mol),
            float(rdMolDescriptors.CalcNumHBD(mol)),
            float(rdMolDescriptors.CalcNumHBA(mol)),
        ])
    return torch.FloatTensor(rows)


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

class PropertyPredictor(nn.Module):
    """
    Multi-output MLP with MC-Dropout uncertainty.

    Architecture:
        2048 → 512 → Dropout(0.3) → 256 → Dropout(0.3) → N_PROPS heads
    """

    def __init__(self, input_dim: int = 2048, dropout: float = 0.3):
        super().__init__()
        self.backbone = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        # One head per property for independent loss scaling
        self.heads = nn.ModuleList([nn.Linear(256, 1) for _ in range(N_PROPS)])

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Returns (batch, N_PROPS) predictions."""
        h = self.backbone(x)
        return torch.cat([head(h) for head in self.heads], dim=1)


# ---------------------------------------------------------------------------
# Normalisation helpers
# ---------------------------------------------------------------------------

@dataclass
class Scaler:
    mean: torch.Tensor
    std:  torch.Tensor

    def transform(self, y: torch.Tensor) -> torch.Tensor:
        return (y - self.mean) / (self.std + 1e-8)

    def inverse(self, y: torch.Tensor) -> torch.Tensor:
        return y * (self.std + 1e-8) + self.mean


def fit_scaler(y: torch.Tensor) -> Scaler:
    return Scaler(mean=y.mean(0), std=y.std(0))


# ---------------------------------------------------------------------------
# MC-Dropout uncertainty
# ---------------------------------------------------------------------------

def mc_property_predict(
    model: PropertyPredictor,
    x: torch.Tensor,
    samples: int = 50,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Run `samples` stochastic forward passes and return (mean, variance).
    Both tensors have shape (batch, N_PROPS).
    """
    model.train()   # keep dropout active
    preds = []
    with torch.no_grad():
        for _ in range(samples):
            preds.append(model(x))
    stack = torch.stack(preds)          # (samples, batch, N_PROPS)
    return stack.mean(0), stack.var(0)


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------

def train(epochs: int = 300, batch_size: int = 256, lr: float = 1e-3) -> PropertyPredictor:
    x_path = os.path.join(DATA_DIR, "X.pt")
    csv_path = os.path.join(DATA_DIR, "molecules_clean.csv")

    if not os.path.exists(x_path) or not os.path.exists(csv_path):
        raise FileNotFoundError(
            "Run src/data/data_loader.py first to generate X.pt and molecules_clean.csv."
        )

    print("Loading features and computing property labels ...")
    X = torch.load(x_path, weights_only=True)           # (N, 2048)
    Y = build_property_labels(csv_path)[:len(X)]        # (N, 6)

    scaler = fit_scaler(Y)
    Y_norm = scaler.transform(Y)

    model     = PropertyPredictor()
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss()

    N = X.shape[0]
    print(f"Training PropertyPredictor: {N} molecules, {epochs} epochs ...")

    for epoch in range(epochs):
        # Mini-batch loop
        perm = torch.randperm(N)
        epoch_loss = 0.0
        for i in range(0, N, batch_size):
            idx = perm[i : i + batch_size]
            xb, yb = X[idx], Y_norm[idx]

            optimizer.zero_grad()
            pred = model(xb)
            loss = criterion(pred, yb)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item() * len(idx)

        if (epoch + 1) % 50 == 0:
            print(f"  Epoch {epoch+1}/{epochs} | MSE: {epoch_loss/N:.6f}")

    # Save model + scaler
    ckpt_path = os.path.join(DATA_DIR, "property_predictor.pt")
    torch.save(
        {"model": model.state_dict(), "scaler_mean": scaler.mean, "scaler_std": scaler.std},
        ckpt_path,
    )
    print(f"\nModel saved → {ckpt_path}")
    return model


# ---------------------------------------------------------------------------
# Convenience: predict for a single SMILES
# ---------------------------------------------------------------------------

def predict_smiles(smiles: str, ckpt_path: str = None) -> dict:
    """
    Load a trained model and return a property dict for a SMILES string.
    Includes mean prediction and uncertainty (std across MC samples).
    """
    from rdkit.Chem import rdFingerprintGenerator

    if ckpt_path is None:
        ckpt_path = os.path.join(DATA_DIR, "property_predictor.pt")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp  = torch.FloatTensor(gen.GetFingerprintAsNumPy(mol)).unsqueeze(0)

    ckpt   = torch.load(ckpt_path, weights_only=True)
    scaler = Scaler(mean=ckpt["scaler_mean"], std=ckpt["scaler_std"])
    model  = PropertyPredictor()
    model.load_state_dict(ckpt["model"])

    mean_norm, var_norm = mc_property_predict(model, fp, samples=50)
    mean = scaler.inverse(mean_norm).squeeze(0)
    std  = (var_norm.sqrt() * (scaler.std + 1e-8)).squeeze(0)

    return {
        name: {"mean": round(mean[i].item(), 4), "std": round(std[i].item(), 4)}
        for i, name in enumerate(PROPERTY_NAMES)
    }


if __name__ == "__main__":
    train()
