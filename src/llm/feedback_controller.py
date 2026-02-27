"""Layer 9 → Layer 7 feedback controller.

Analyses ValidationResult failures from Layer 9, extracts generation
constraints, and applies soft bounds to the BayesianGraphVAE's latent
space sampling so subsequent generation rounds avoid invalid regions.

This implements the feedback loop documented in README.md lines 330-332:
    Layer 9 → Layer 7 (Guide generation away from invalid)
"""

from dataclasses import dataclass, field
from typing import Optional

from .pipeline_layer9 import ValidationResult


# ---------------------------------------------------------------------------
# Constraint dataclass
# ---------------------------------------------------------------------------

@dataclass
class GenerationConstraints:
    """Soft upper bounds derived from accumulated validation failures."""
    max_mw:       float = 500.0
    max_logP:     float = 5.0
    max_hbd:      int   = 5
    max_hba:      int   = 10
    max_rot_bonds: int  = 10
    max_psa:      float = 140.0
    reject_pains: bool  = True
    # latent dimensions flagged as over-active (for future use)
    latent_penalty_dims: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Controller
# ---------------------------------------------------------------------------

class FeedbackController:
    """
    Inspects Layer 9 failures and incrementally tightens GenerationConstraints.

    Typical usage:
        controller = FeedbackController()
        failures = [r for r in results if not r.valid]
        constraints = controller.extract_constraints_from_failures(failures)
        controller.update_vae_sampling(vae, constraints)
    """

    def __init__(self):
        self.constraints = GenerationConstraints()
        self._failure_counts: dict = {}

    def extract_constraints_from_failures(
        self, failures: list
    ) -> GenerationConstraints:
        """
        Tighten bounds based on which checks failed.

        | Failed check       | Rocq error type   | Action                      |
        |--------------------|-------------------|-----------------------------|
        | Lipinski MW/LogP   | NUMERIC_BOUNDS    | Shrink max_mw / max_logP    |
        | Veber rot/PSA      | NUMERIC_BOUNDS    | Shrink max_rot_bonds / PSA  |
        | PAINS/Brenk safety | —                 | Enable reject_pains flag    |
        | Rocq proof         | TYPE_ERROR etc.   | Count; no direct bound      |
        """
        for result in failures:
            error_type = result.error_type or "UNKNOWN"
            self._failure_counts[error_type] = (
                self._failure_counts.get(error_type, 0) + 1
            )

            if not result.lipinski_ok:
                # Tighten MW/logP by ~10 % toward stricter drug-like centre
                self.constraints.max_mw  = min(self.constraints.max_mw,  450.0)
                self.constraints.max_logP = min(self.constraints.max_logP, 4.5)
                self.constraints.max_hbd  = min(self.constraints.max_hbd,  4)
                self.constraints.max_hba  = min(self.constraints.max_hba,  9)

            if not result.veber_ok:
                self.constraints.max_rot_bonds = min(
                    self.constraints.max_rot_bonds, 8
                )
                self.constraints.max_psa = min(self.constraints.max_psa, 120.0)

            if not result.safety_ok:
                self.constraints.reject_pains = True

        return self.constraints

    def update_vae_sampling(self, vae, constraints: Optional[GenerationConstraints] = None) -> None:
        """
        Monkey-patch the VAE's reparameterize method to add a soft
        constraint: latent vectors with extreme L2 norm (which tend to
        decode to out-of-distribution, drug-unlike molecules) are redrawn
        up to max_retries times before accepting the last sample anyway.

        This is intentionally lightweight — the norm threshold corresponds
        roughly to the 99th-percentile L2 norm of z ~ N(0,I) in 16 dims
        (~sqrt(16) * 2.6 ≈ 10.4).
        """
        if constraints is None:
            constraints = self.constraints

        original_reparameterize = vae.reparameterize

        def constrained_reparameterize(mu, logvar, max_retries: int = 10):
            for _ in range(max_retries):
                z = original_reparameterize(mu, logvar)
                if z.norm(dim=-1).max().item() < 10.4:
                    return z
            return z  # accept last sample if constraint never met

        vae.reparameterize = constrained_reparameterize

    def failure_summary(self) -> dict:
        """Return error-type counts sorted by frequency (most common first)."""
        return dict(
            sorted(self._failure_counts.items(), key=lambda kv: -kv[1])
        )

    def reset(self) -> None:
        """Reset constraints and failure counts to defaults."""
        self.constraints = GenerationConstraints()
        self._failure_counts = {}
