"""Molecule generation benchmark suite.

Provides the following benchmark categories:
  • VUN     – Validity / Uniqueness / Novelty
  • Dist    – Distribution alignment (KL, Wasserstein, MMD) for QED, LogP, MW, etc.
  • Div     – Diversity & mode-collapse (Tanimoto, scaffold diversity, duplicates)
  • FCD     – Fréchet ChemNet Distance
  • Scaffold – Scaffold-split generalization on unseen Bemis–Murcko scaffolds
"""

from .molecule_benchmarks import (
    MoleculeBenchmarks,
    BenchmarkResult,
    compute_vun,
    compute_distribution_alignment,
    compute_diversity,
    compute_fcd,
    compute_scaffold_generalization,
)

__all__ = [
    "MoleculeBenchmarks",
    "BenchmarkResult",
    "compute_vun",
    "compute_distribution_alignment",
    "compute_diversity",
    "compute_fcd",
    "compute_scaffold_generalization",
]
