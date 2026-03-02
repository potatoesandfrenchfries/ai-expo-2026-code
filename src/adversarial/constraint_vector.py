"""Constraint violation vector — maps structured Rocq bridge output to a
differentiable-compatible penalty signal for RL reward shaping.

Public API:
    ConstraintViolationVector             — per-constraint violation state
    ConstraintViolationVector.from_bridge_result(result) → instance
    ConstraintViolationVector.penalty_score(weights)     → float
    load_weights(config_path)             → dict
    default_weights()                     → dict

The vector has one slot per constraint type:
    mw, logp, hbd, hba, rot, psa, pains, sa_score, toxicity

Each slot carries:
    violated  — bool (1 if the constraint failed, 0 otherwise)
    magnitude — float (how far outside the bound; 0 for satisfied constraints)

The scalar penalty score is:
    penalty = sum_i( weight_i * violated_i * clip(magnitude_i, max_clip) )

PAINS and toxicity carry higher default weights than marginal Lipinski
violations (see config/constraint_weights.yaml).
"""

import os
from dataclasses import dataclass, field
from typing import Optional

# ---------------------------------------------------------------------------
# Constraint slot names (canonical order)
# ---------------------------------------------------------------------------

CONSTRAINT_TYPES = [
    "mw",
    "logp",
    "hbd",
    "hba",
    "rot",
    "psa",
    "pains",
    "sa_score",
    "toxicity",
]

_DEFAULT_WEIGHTS = {
    "mw":       1.0,
    "logp":     1.5,
    "hbd":      1.0,
    "hba":      1.0,
    "rot":      0.8,
    "psa":      0.8,
    "pains":    5.0,
    "sa_score": 2.0,
    "toxicity": 5.0,
}

_DEFAULT_MAGNITUDE_CLIP = 10.0


# ---------------------------------------------------------------------------
# ConstraintViolationVector
# ---------------------------------------------------------------------------

@dataclass
class ConstraintViolationVector:
    """Per-constraint violation state for a single molecule.

    Attributes
    ----------
    violated : dict[str, bool]
        1.0 if the constraint failed, 0.0 if satisfied.
    magnitude : dict[str, float]
        How far the property is outside its bound (0.0 for satisfied).
    all_passed : bool
        True only if every constraint slot is satisfied.
    """

    violated:   dict = field(default_factory=lambda: {c: False for c in CONSTRAINT_TYPES})
    magnitude:  dict = field(default_factory=lambda: {c: 0.0   for c in CONSTRAINT_TYPES})
    all_passed: bool = True

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_bridge_result(cls, result: dict) -> "ConstraintViolationVector":
        """Build a vector from a rocq_bridge.check_molecule() result dict."""
        vec = cls()
        vec.all_passed = result.get("passed", True)

        ctype = result.get("constraint_type")
        mag   = result.get("violation_magnitude") or 0.0

        if ctype and ctype in CONSTRAINT_TYPES:
            vec.violated[ctype]  = True
            vec.magnitude[ctype] = float(mag)

        return vec

    @classmethod
    def from_descriptor_dict(cls, desc: dict) -> "ConstraintViolationVector":
        """Build a vector directly from a descriptor dict (no Rocq needed).

        Useful for training data labelling when Rocq is unavailable.
        `desc` must have keys: mw, logp, hbd, hba, rot, psa.
        """
        vec = cls()
        checks = [
            ("mw",   max(0.0, desc.get("mw",   0) - 500)),
            ("logp", max(0.0, desc.get("logp",  0) - 5)),
            ("hbd",  max(0.0, float(desc.get("hbd", 0) - 5))),
            ("hba",  max(0.0, float(desc.get("hba", 0) - 10))),
            ("rot",  max(0.0, float(desc.get("rot", 0) - 10))),
            ("psa",  max(0.0, desc.get("psa",  0) - 140)),
        ]
        for ctype, mag in checks:
            if mag > 0:
                vec.violated[ctype]  = True
                vec.magnitude[ctype] = mag
                vec.all_passed = False

        # Binary flags (caller can set via mutate after construction)
        for bin_type in ("pains", "sa_score", "toxicity"):
            if desc.get(bin_type, False):
                vec.violated[bin_type]  = True
                vec.magnitude[bin_type] = 1.0
                vec.all_passed = False

        return vec

    # ------------------------------------------------------------------
    # Penalty computation
    # ------------------------------------------------------------------

    def penalty_score(
        self,
        weights: Optional[dict] = None,
        magnitude_clip: float = _DEFAULT_MAGNITUDE_CLIP,
    ) -> float:
        """Compute the weighted scalar penalty.

        penalty = sum_i( w_i * violated_i * clip(magnitude_i, magnitude_clip) )

        Parameters
        ----------
        weights : dict, optional
            Per-constraint weights; defaults to _DEFAULT_WEIGHTS.
        magnitude_clip : float
            Maximum magnitude contribution per constraint.
        """
        w = weights if weights is not None else _DEFAULT_WEIGHTS
        total = 0.0
        for ctype in CONSTRAINT_TYPES:
            if self.violated[ctype]:
                mag = min(float(self.magnitude[ctype]), magnitude_clip)
                total += w.get(ctype, 1.0) * mag
        return total

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """Flat dict for logging / metrics."""
        out = {"all_passed": self.all_passed}
        for ctype in CONSTRAINT_TYPES:
            out[f"violated_{ctype}"]  = self.violated[ctype]
            out[f"magnitude_{ctype}"] = self.magnitude[ctype]
        return out

    def violated_list(self) -> list:
        """Return a list of constraint types that are violated."""
        return [c for c in CONSTRAINT_TYPES if self.violated[c]]

    def as_binary_vector(self) -> list:
        """Return a list of 0/1 integers in CONSTRAINT_TYPES order."""
        return [int(self.violated[c]) for c in CONSTRAINT_TYPES]

    def as_magnitude_vector(self) -> list:
        """Return a list of violation magnitudes in CONSTRAINT_TYPES order."""
        return [self.magnitude[c] for c in CONSTRAINT_TYPES]

    def __repr__(self) -> str:
        vs = [c for c in CONSTRAINT_TYPES if self.violated[c]]
        status = "PASS" if self.all_passed else f"FAIL({', '.join(vs)})"
        return f"ConstraintViolationVector({status})"


# ---------------------------------------------------------------------------
# Weight loading
# ---------------------------------------------------------------------------

def load_weights(config_path: Optional[str] = None) -> dict:
    """Load constraint weights from a YAML config file.

    Falls back to default weights if the file is absent or yaml is
    not installed.
    """
    if config_path is None:
        repo_root = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..")
        )
        config_path = os.path.join(repo_root, "config", "constraint_weights.yaml")

    try:
        import yaml
        with open(config_path) as f:
            cfg = yaml.safe_load(f)
        weights = {k: float(v) for k, v in cfg.get("weights", {}).items()}
        # Fill in any missing keys with defaults
        for k, v in _DEFAULT_WEIGHTS.items():
            weights.setdefault(k, v)
        return weights
    except Exception:
        return dict(_DEFAULT_WEIGHTS)


def default_weights() -> dict:
    """Return the built-in default weights (no file I/O)."""
    return dict(_DEFAULT_WEIGHTS)
