"""Molecule generation benchmark suite.

Five benchmark categories mirroring the current state of the art for evaluating
generative molecule models:

1. VUN – Validity / Uniqueness / Novelty
   • % chemically valid SMILES
   • % unique canonical SMILES among the valid ones
   • novelty against train set and holdout set separately

2. Distribution alignment
   • Per-property KL divergence and Wasserstein-1 distance (QED, LogP, MW, TPSA,
     SA score, ring count, heteroatom count)
   • Maximum Mean Discrepancy (MMD) with RBF kernel over the joint property vector

3. Diversity & mode collapse
   • Internal mean pairwise Tanimoto distance (Morgan FP)
   • Scaffold diversity: fraction of unique Bemis–Murcko scaffolds
   • % duplicate scaffolds

4. Fréchet ChemNet Distance (FCD)
   • Computes activation statistics from the ChemNet hidden layer and evaluates
     the Fréchet distance between generated and reference sets.
   • Falls back to a lightweight RDKit-fingerprint-based Fréchet surrogate when
     the ``fcd`` package is not installed.

5. Scaffold-split generalization
   • Splits a reference corpus by Bemis–Murcko scaffold into seen / unseen halves
   • Reports novelty and validity/uniqueness metrics for generated molecules that
     map to the *unseen* scaffold partition

Usage
─────
    from src.benchmarks import MoleculeBenchmarks

    benchmarks = MoleculeBenchmarks(
        train_smiles=train_list,
        holdout_smiles=holdout_list,
    )
    result = benchmarks.evaluate(generated_smiles)
    result.print_summary()
    d = result.to_dict()
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import (
        Descriptors,
        QED,
        rdMolDescriptors,
    )
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False
    warnings.warn("RDKit is not installed; most benchmarks will be unavailable.", stacklevel=2)

try:
    from rdkit.Chem import RDConfig
    import sys
    import os

    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
    import sascorer  # type: ignore

    _SA_AVAILABLE = True
except Exception:
    _SA_AVAILABLE = False

try:
    import fcd as _fcd_module  # type: ignore

    _FCD_AVAILABLE = True
except ImportError:
    _FCD_AVAILABLE = False

try:
    from scipy import stats as _scipy_stats  # type: ignore

    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Helper: canonicalise and validate SMILES
# ---------------------------------------------------------------------------

def _canonical(smiles: str) -> Optional[str]:
    """Return canonical SMILES or None if chemically invalid."""
    if not _RDKIT_AVAILABLE:
        return smiles  # cannot validate without RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def _validate_batch(smiles_list: Sequence[str]) -> Tuple[List[str], List[str]]:
    """Return (valid_canonical_list, invalid_list)."""
    valid, invalid = [], []
    for s in smiles_list:
        c = _canonical(s)
        if c is not None:
            valid.append(c)
        else:
            invalid.append(s)
    return valid, invalid


# ---------------------------------------------------------------------------
# Helper: Morgan fingerprints
# ---------------------------------------------------------------------------

def _morgan_fp(smiles: str, radius: int = 2, n_bits: int = 2048):
    """Return an RDKit ExplicitBitVect Morgan fingerprint, or None."""
    if not _RDKIT_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def _fp_matrix(smiles_list: Sequence[str], radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    """Return (N, n_bits) binary numpy array of Morgan fingerprints."""
    rows = []
    for s in smiles_list:
        fp = _morgan_fp(s, radius=radius, n_bits=n_bits)
        if fp is not None:
            arr = np.zeros(n_bits, dtype=np.float32)
            DataStructs.ConvertToNumpyArray(fp, arr)
            rows.append(arr)
    if not rows:
        return np.zeros((0, n_bits), dtype=np.float32)
    return np.stack(rows, axis=0)


# ---------------------------------------------------------------------------
# Helper: per-molecule property computation
# ---------------------------------------------------------------------------

def _compute_properties(smiles: str) -> Optional[Dict[str, float]]:
    """Compute QED, LogP, MW, TPSA, SA, ring_count, heteroatom_count."""
    if not _RDKIT_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    sa = float("nan")
    if _SA_AVAILABLE:
        try:
            sa = sascorer.calculateScore(mol)
        except Exception:
            pass

    return {
        "qed":            QED.qed(mol),
        "logp":           Descriptors.MolLogP(mol),
        "mw":             Descriptors.MolWt(mol),
        "tpsa":           rdMolDescriptors.CalcTPSA(mol),
        "sa":             sa,
        "ring_count":     rdMolDescriptors.CalcNumRings(mol),
        "heteroatom_count": sum(
            1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6)
        ),
    }


def _property_arrays(smiles_list: Sequence[str]) -> Dict[str, np.ndarray]:
    """Return dict of property → 1-D float array, skipping failed computations."""
    records = [_compute_properties(s) for s in smiles_list]
    records = [r for r in records if r is not None]
    if not records:
        return {}
    keys = list(records[0].keys())
    out: Dict[str, np.ndarray] = {}
    for k in keys:
        vals = np.array([r[k] for r in records], dtype=np.float64)
        finite = vals[np.isfinite(vals)]
        if len(finite) > 0:
            out[k] = finite
    return out


# ---------------------------------------------------------------------------
# Helper: Bemis–Murcko scaffold
# ---------------------------------------------------------------------------

def _bemis_murcko(smiles: str) -> Optional[str]:
    """Return canonical Bemis–Murcko scaffold SMILES, or None."""
    if not _RDKIT_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold, canonical=True)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Distribution-alignment helpers
# ---------------------------------------------------------------------------

def _kl_divergence_binned(p: np.ndarray, q: np.ndarray, bins: int = 50) -> float:
    """KL(P ‖ Q) estimated from samples via equal-width histogram binning."""
    lo = min(p.min(), q.min())
    hi = max(p.max(), q.max())
    if lo == hi:
        return 0.0
    edges = np.linspace(lo, hi, bins + 1)
    p_hist, _ = np.histogram(p, bins=edges, density=True)
    q_hist, _ = np.histogram(q, bins=edges, density=True)
    # Smooth to avoid log(0)
    eps = 1e-10
    p_hist = p_hist + eps
    q_hist = q_hist + eps
    p_hist /= p_hist.sum()
    q_hist /= q_hist.sum()
    return float(np.sum(p_hist * np.log(p_hist / q_hist)))


def _wasserstein1(p: np.ndarray, q: np.ndarray) -> float:
    """Wasserstein-1 (Earth-Mover) distance between two 1-D sample arrays."""
    if _SCIPY_AVAILABLE:
        return float(_scipy_stats.wasserstein_distance(p, q))
    # Fallback: sort-based exact computation
    p_sorted = np.sort(p)
    q_sorted = np.sort(q)
    # Interpolate to common length
    n = max(len(p_sorted), len(q_sorted))
    p_i = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(p_sorted)), p_sorted)
    q_i = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(q_sorted)), q_sorted)
    return float(np.mean(np.abs(p_i - q_i)))


def _mmd_rbf(X: np.ndarray, Y: np.ndarray, gamma: float = 1.0) -> float:
    """Maximum Mean Discrepancy with RBF kernel, O(N²) implementation."""
    if X.shape[0] == 0 or Y.shape[0] == 0:
        return float("nan")

    def rbf(A: np.ndarray, B: np.ndarray) -> float:
        sq = (
            np.sum(A ** 2, axis=1, keepdims=True)
            + np.sum(B ** 2, axis=1, keepdims=True).T
            - 2.0 * A @ B.T
        )
        return float(np.mean(np.exp(-gamma * sq)))

    # Subsample for efficiency if sets are large
    max_n = 500
    if X.shape[0] > max_n:
        idx = np.random.choice(X.shape[0], max_n, replace=False)
        X = X[idx]
    if Y.shape[0] > max_n:
        idx = np.random.choice(Y.shape[0], max_n, replace=False)
        Y = Y[idx]

    return rbf(X, X) - 2.0 * rbf(X, Y) + rbf(Y, Y)


# ---------------------------------------------------------------------------
# FCD helpers
# ---------------------------------------------------------------------------

def _fcd_score(generated: Sequence[str], reference: Sequence[str]) -> float:
    """Compute Fréchet ChemNet Distance.

    Uses the ``fcd`` package (https://github.com/insightsengineering/fcd) when
    available; otherwise falls back to a surrogate computed from Morgan
    fingerprint activation statistics.
    """
    if _FCD_AVAILABLE:
        try:
            return float(_fcd_module.get_fcd(list(generated), list(reference)))
        except Exception as exc:
            warnings.warn(f"fcd package raised {exc}; using fingerprint surrogate.", stacklevel=3)

    # Fingerprint-based Fréchet surrogate
    gen_fps = _fp_matrix(list(generated))
    ref_fps = _fp_matrix(list(reference))
    if gen_fps.shape[0] < 2 or ref_fps.shape[0] < 2:
        return float("nan")

    mu_g, sigma_g = gen_fps.mean(axis=0), np.cov(gen_fps, rowvar=False)
    mu_r, sigma_r = ref_fps.mean(axis=0), np.cov(ref_fps, rowvar=False)

    diff = mu_g - mu_r
    # Compute sqrt of matrix product via eigendecomposition (more stable than sqrtm)
    try:
        # sigma_g @ sigma_r
        product = sigma_g @ sigma_r
        eigvals = np.linalg.eigvalsh(product)
        sqrt_trace = float(np.sum(np.sqrt(np.maximum(eigvals, 0.0))))
    except np.linalg.LinAlgError:
        sqrt_trace = 0.0

    fid = float(
        diff @ diff
        + np.trace(sigma_g)
        + np.trace(sigma_r)
        - 2.0 * sqrt_trace
    )
    return max(fid, 0.0)


# ---------------------------------------------------------------------------
# VUN
# ---------------------------------------------------------------------------

def compute_vun(
    generated_smiles: Sequence[str],
    train_smiles: Sequence[str],
    holdout_smiles: Optional[Sequence[str]] = None,
) -> Dict[str, float]:
    """Compute Validity / Uniqueness / Novelty metrics.

    Parameters
    ----------
    generated_smiles : list of str
        Raw SMILES produced by the generative model.
    train_smiles : list of str
        Canonical SMILES used during training.
    holdout_smiles : list of str, optional
        Canonical SMILES from a held-out reference set.

    Returns
    -------
    dict with keys:
        validity          – fraction of chemically valid SMILES
        uniqueness        – fraction of unique among valid
        novelty_vs_train  – fraction of unique valid SMILES not in train set
        novelty_vs_holdout – fraction not in holdout (NaN if not provided)
    """
    n_total = len(generated_smiles)
    if n_total == 0:
        return {
            "validity": float("nan"),
            "uniqueness": float("nan"),
            "novelty_vs_train": float("nan"),
            "novelty_vs_holdout": float("nan"),
        }

    valid_list, _ = _validate_batch(generated_smiles)
    validity = len(valid_list) / n_total

    if not valid_list:
        return {
            "validity": validity,
            "uniqueness": float("nan"),
            "novelty_vs_train": float("nan"),
            "novelty_vs_holdout": float("nan"),
        }

    unique_set = set(valid_list)
    uniqueness = len(unique_set) / len(valid_list)

    train_set = set(train_smiles)
    novelty_train = sum(1 for s in unique_set if s not in train_set) / len(unique_set)

    novelty_holdout = float("nan")
    if holdout_smiles is not None:
        holdout_set = set(holdout_smiles)
        novelty_holdout = sum(1 for s in unique_set if s not in holdout_set) / len(unique_set)

    return {
        "validity": validity,
        "uniqueness": uniqueness,
        "novelty_vs_train": novelty_train,
        "novelty_vs_holdout": novelty_holdout,
    }


# ---------------------------------------------------------------------------
# Distribution alignment
# ---------------------------------------------------------------------------

def compute_distribution_alignment(
    generated_smiles: Sequence[str],
    reference_smiles: Sequence[str],
    properties: Optional[List[str]] = None,
    mmd_gamma: float = 1.0,
) -> Dict[str, float]:
    """Compute per-property KL divergence, Wasserstein-1, and MMD.

    Parameters
    ----------
    generated_smiles : list of str
    reference_smiles : list of str
    properties : list of str, optional
        Subset of ['qed', 'logp', 'mw', 'tpsa', 'sa', 'ring_count', 'heteroatom_count'].
        Defaults to all.
    mmd_gamma : float
        RBF kernel bandwidth for MMD (on the joint fingerprint space).

    Returns
    -------
    dict with keys like:
        kl_qed, kl_logp, …
        wasserstein_qed, wasserstein_logp, …
        mmd_fingerprint
    """
    _all_props = ["qed", "logp", "mw", "tpsa", "sa", "ring_count", "heteroatom_count"]
    if properties is None:
        properties = _all_props

    gen_props = _property_arrays(generated_smiles)
    ref_props = _property_arrays(reference_smiles)

    result: Dict[str, float] = {}
    for prop in properties:
        if prop not in gen_props or prop not in ref_props:
            result[f"kl_{prop}"] = float("nan")
            result[f"wasserstein_{prop}"] = float("nan")
            continue
        g = gen_props[prop]
        r = ref_props[prop]
        result[f"kl_{prop}"] = _kl_divergence_binned(g, r)
        result[f"wasserstein_{prop}"] = _wasserstein1(g, r)

    # MMD on fingerprint space
    gen_fps = _fp_matrix(list(generated_smiles))
    ref_fps = _fp_matrix(list(reference_smiles))
    result["mmd_fingerprint"] = _mmd_rbf(gen_fps, ref_fps, gamma=mmd_gamma)

    return result


# ---------------------------------------------------------------------------
# Diversity & mode collapse
# ---------------------------------------------------------------------------

def compute_diversity(generated_smiles: Sequence[str]) -> Dict[str, float]:
    """Compute internal diversity metrics for a set of generated molecules.

    Parameters
    ----------
    generated_smiles : list of str
        Canonicalised valid SMILES (the caller should pre-filter invalid ones).

    Returns
    -------
    dict with keys:
        tanimoto_diversity   – mean pairwise 1 − Tanimoto similarity
        scaffold_diversity   – fraction of unique Bemis–Murcko scaffolds
        duplicate_scaffold_frac – fraction of molecules sharing a scaffold
    """
    if not _RDKIT_AVAILABLE:
        return {
            "tanimoto_diversity": float("nan"),
            "scaffold_diversity": float("nan"),
            "duplicate_scaffold_frac": float("nan"),
        }

    smiles = list(generated_smiles)
    n = len(smiles)
    if n == 0:
        return {
            "tanimoto_diversity": float("nan"),
            "scaffold_diversity": float("nan"),
            "duplicate_scaffold_frac": float("nan"),
        }

    # ---- Tanimoto diversity (subsample for large sets) ----
    fps = [_morgan_fp(s) for s in smiles]
    fps = [fp for fp in fps if fp is not None]

    max_pairs = 5000
    if len(fps) > 1:
        # Randomly sample pairs
        rng = np.random.default_rng(42)
        idx_pairs = rng.integers(0, len(fps), size=(min(max_pairs, len(fps) * (len(fps) - 1) // 2), 2))
        # Ensure i != j
        mask = idx_pairs[:, 0] != idx_pairs[:, 1]
        idx_pairs = idx_pairs[mask]
        if len(idx_pairs) == 0:
            tanimoto_div = float("nan")
        else:
            sims = [
                DataStructs.TanimotoSimilarity(fps[i], fps[j])
                for i, j in idx_pairs
            ]
            tanimoto_div = 1.0 - float(np.mean(sims))
    else:
        tanimoto_div = float("nan")

    # ---- Scaffold diversity ----
    scaffolds = [_bemis_murcko(s) for s in smiles]
    scaffolds = [sc for sc in scaffolds if sc is not None]
    if scaffolds:
        unique_scaffolds = set(scaffolds)
        scaffold_diversity = len(unique_scaffolds) / len(scaffolds)

        # % molecules whose scaffold appears more than once (duplicate)
        from collections import Counter
        scaffold_counts = Counter(scaffolds)
        dup_count = sum(1 for sc in scaffolds if scaffold_counts[sc] > 1)
        duplicate_scaffold_frac = dup_count / len(scaffolds)
    else:
        scaffold_diversity = float("nan")
        duplicate_scaffold_frac = float("nan")

    return {
        "tanimoto_diversity": tanimoto_div,
        "scaffold_diversity": scaffold_diversity,
        "duplicate_scaffold_frac": duplicate_scaffold_frac,
    }


# ---------------------------------------------------------------------------
# FCD
# ---------------------------------------------------------------------------

def compute_fcd(
    generated_smiles: Sequence[str],
    reference_smiles: Sequence[str],
) -> Dict[str, float]:
    """Compute Fréchet ChemNet Distance between generated and reference sets.

    Returns
    -------
    dict with key:
        fcd  – scalar distance (lower = better, 0 = identical distributions)
    """
    return {"fcd": _fcd_score(generated_smiles, reference_smiles)}


# ---------------------------------------------------------------------------
# Scaffold-split generalisation
# ---------------------------------------------------------------------------

def compute_scaffold_generalization(
    generated_smiles: Sequence[str],
    reference_smiles: Sequence[str],
    unseen_fraction: float = 0.5,
    seed: int = 42,
) -> Dict[str, float]:
    """Evaluate novelty and quality on molecules with unseen Bemis–Murcko scaffolds.

    The reference corpus is split by scaffold into a *seen* half and an *unseen*
    half (by number of scaffolds, not molecules).  Generated molecules whose
    canonical scaffold falls in the unseen partition are then evaluated for
    validity and novelty.

    Parameters
    ----------
    generated_smiles : list of str
    reference_smiles : list of str
        Full reference corpus used to define seen/unseen scaffold split.
    unseen_fraction : float
        Fraction of unique reference scaffolds to designate as unseen (default 0.5).
    seed : int
        Random seed for scaffold assignment.

    Returns
    -------
    dict with keys:
        unseen_scaffold_count      – number of unique unseen scaffolds
        generated_on_unseen_count  – how many generated SMILES map to unseen scaffolds
        novelty_on_unseen          – fraction of those not in reference corpus
        validity_on_unseen         – fraction that are chemically valid
        uniqueness_on_unseen       – fraction that are unique among those matched
    """
    if not _RDKIT_AVAILABLE:
        return {
            "unseen_scaffold_count": float("nan"),
            "generated_on_unseen_count": float("nan"),
            "novelty_on_unseen": float("nan"),
            "validity_on_unseen": float("nan"),
            "uniqueness_on_unseen": float("nan"),
        }

    # Build scaffold → [smiles] map for the reference corpus
    ref_scaffold_map: Dict[str, List[str]] = {}
    for s in reference_smiles:
        sc = _bemis_murcko(s)
        if sc is None:
            sc = "__no_scaffold__"
        ref_scaffold_map.setdefault(sc, []).append(s)

    all_scaffolds = list(ref_scaffold_map.keys())
    rng = np.random.default_rng(seed)
    rng.shuffle(all_scaffolds)
    n_unseen = max(1, int(len(all_scaffolds) * unseen_fraction))
    unseen_scaffolds: set = set(all_scaffolds[:n_unseen])
    ref_canonical: set = set(
        c
        for s in reference_smiles
        for c in [_canonical(s)]
        if c is not None
    )

    # Filter generated molecules to those mapping to unseen scaffolds
    on_unseen_raw: List[str] = []
    on_unseen_valid: List[str] = []
    for s in generated_smiles:
        sc = _bemis_murcko(s)
        if sc is None:
            sc = "__no_scaffold__"
        if sc in unseen_scaffolds:
            on_unseen_raw.append(s)
            c = _canonical(s)
            if c is not None:
                on_unseen_valid.append(c)

    n_raw = len(on_unseen_raw)
    if n_raw == 0:
        return {
            "unseen_scaffold_count": len(unseen_scaffolds),
            "generated_on_unseen_count": 0,
            "novelty_on_unseen": float("nan"),
            "validity_on_unseen": float("nan"),
            "uniqueness_on_unseen": float("nan"),
        }

    validity_on_unseen = len(on_unseen_valid) / n_raw
    unique_on_unseen = set(on_unseen_valid)
    uniqueness_on_unseen = len(unique_on_unseen) / len(on_unseen_valid) if on_unseen_valid else float("nan")
    novelty_on_unseen = (
        sum(1 for s in unique_on_unseen if s not in ref_canonical) / len(unique_on_unseen)
        if unique_on_unseen else float("nan")
    )

    return {
        "unseen_scaffold_count": len(unseen_scaffolds),
        "generated_on_unseen_count": n_raw,
        "novelty_on_unseen": novelty_on_unseen,
        "validity_on_unseen": validity_on_unseen,
        "uniqueness_on_unseen": uniqueness_on_unseen,
    }


# ---------------------------------------------------------------------------
# BenchmarkResult
# ---------------------------------------------------------------------------

@dataclass
class BenchmarkResult:
    """Aggregated results from all benchmark categories."""

    vun: Dict[str, float] = field(default_factory=dict)
    distribution: Dict[str, float] = field(default_factory=dict)
    diversity: Dict[str, float] = field(default_factory=dict)
    fcd: Dict[str, float] = field(default_factory=dict)
    scaffold_generalization: Dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, float]:
        """Flatten all metrics into a single dict."""
        out: Dict[str, float] = {}
        for prefix, d in [
            ("vun", self.vun),
            ("dist", self.distribution),
            ("div", self.diversity),
            ("fcd", self.fcd),
            ("scaffold", self.scaffold_generalization),
        ]:
            for k, v in d.items():
                out[f"{prefix}/{k}"] = v
        return out

    def print_summary(self) -> None:
        """Print a human-readable summary to stdout."""
        sections = [
            ("VUN", self.vun),
            ("Distribution alignment", self.distribution),
            ("Diversity", self.diversity),
            ("FCD", self.fcd),
            ("Scaffold generalisation", self.scaffold_generalization),
        ]
        for title, metrics in sections:
            print(f"\n{'─' * 50}")
            print(f"  {title}")
            print(f"{'─' * 50}")
            for k, v in metrics.items():
                if isinstance(v, float) and not math.isnan(v):
                    print(f"  {k:<35s} {v:.4f}")
                elif isinstance(v, float):
                    print(f"  {k:<35s} NaN")
                else:
                    print(f"  {k:<35s} {v}")


# ---------------------------------------------------------------------------
# MoleculeBenchmarks orchestrator
# ---------------------------------------------------------------------------

class MoleculeBenchmarks:
    """Run all benchmark categories against a fixed train/holdout reference.

    Parameters
    ----------
    train_smiles : list of str
        SMILES used to train the generative model.
    holdout_smiles : list of str, optional
        A held-out reference split (e.g. from ChEMBL or ZINC test partition).
    mmd_gamma : float
        Bandwidth for the MMD RBF kernel.
    scaffold_unseen_fraction : float
        Fraction of scaffolds designated as "unseen" for the scaffold-split test.
    scaffold_seed : int
        RNG seed for scaffold partitioning.
    """

    def __init__(
        self,
        train_smiles: Sequence[str],
        holdout_smiles: Optional[Sequence[str]] = None,
        mmd_gamma: float = 1.0,
        scaffold_unseen_fraction: float = 0.5,
        scaffold_seed: int = 42,
    ):
        self._train = list(train_smiles)
        self._holdout = list(holdout_smiles) if holdout_smiles is not None else None
        self._mmd_gamma = mmd_gamma
        self._scaffold_unseen_fraction = scaffold_unseen_fraction
        self._scaffold_seed = scaffold_seed

        # Pre-canonicalise reference sets for novelty checks
        self._train_canonical = [c for s in self._train for c in [_canonical(s)] if c]
        if self._holdout is not None:
            self._holdout_canonical = [c for s in self._holdout for c in [_canonical(s)] if c]
        else:
            self._holdout_canonical = None

        # Reference set for distribution/scaffold benchmarks: train + holdout combined
        self._reference = self._train_canonical + (self._holdout_canonical or [])

    def evaluate(
        self,
        generated_smiles: Sequence[str],
        skip: Optional[List[str]] = None,
    ) -> BenchmarkResult:
        """Run all benchmarks on ``generated_smiles``.

        Parameters
        ----------
        generated_smiles : list of str
        skip : list of str, optional
            Category names to skip: 'vun', 'distribution', 'diversity', 'fcd',
            'scaffold'.

        Returns
        -------
        BenchmarkResult
        """
        skip = set(skip or [])

        # Pre-filter to valid canonical SMILES for most sub-benchmarks
        valid_list, _ = _validate_batch(generated_smiles)

        vun_metrics: Dict[str, float] = {}
        if "vun" not in skip:
            vun_metrics = compute_vun(
                generated_smiles,
                self._train_canonical,
                self._holdout_canonical,
            )

        dist_metrics: Dict[str, float] = {}
        if "distribution" not in skip and self._reference and valid_list:
            dist_metrics = compute_distribution_alignment(
                valid_list,
                self._reference,
                mmd_gamma=self._mmd_gamma,
            )

        div_metrics: Dict[str, float] = {}
        if "diversity" not in skip and valid_list:
            div_metrics = compute_diversity(valid_list)

        fcd_metrics: Dict[str, float] = {}
        if "fcd" not in skip and self._reference and valid_list:
            fcd_metrics = compute_fcd(valid_list, self._reference)

        scaffold_metrics: Dict[str, float] = {}
        if "scaffold" not in skip and self._reference:
            scaffold_metrics = compute_scaffold_generalization(
                list(generated_smiles),
                self._reference,
                unseen_fraction=self._scaffold_unseen_fraction,
                seed=self._scaffold_seed,
            )

        return BenchmarkResult(
            vun=vun_metrics,
            distribution=dist_metrics,
            diversity=div_metrics,
            fcd=fcd_metrics,
            scaffold_generalization=scaffold_metrics,
        )
