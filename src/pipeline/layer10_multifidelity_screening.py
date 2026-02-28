"""Layer 10: Multi-fidelity molecular screening.

Ranks Layer 9 validated candidates across three increasing fidelity
levels using RDKit-only metrics (no MD/QM infrastructure in Phase 1).
Each level acts as a progressively tighter filter:

  Level 1 (ML gate)        : Layer 9 confidence score [0, 1]
  Level 2 (physicochemical): QED + LogP sweet-spot [1–3] + Bertz complexity
  Level 3 (tractability)   : MW sweet-spot [300–450 Da] + ring count

Final composite score: 0.30 × L1 + 0.40 × L2 + 0.30 × L3

Future integration points (stubbed with TODO markers):
  Level 2 → AutoDock Vina binding energy (replaces physicochemical score)
  Level 3 → GROMACS/OpenMM MM-GBSA ΔG_bind (replaces tractability proxy)

Data contract:
  Input : list[ValidationResult]  from Layer 9 (pass AND fail entries)
  Output: Layer10Result

Usage:
    from src.pipeline.layer10_multifidelity_screening import Layer10MultiFidelityScreening
    result = Layer10MultiFidelityScreening().run(validated_candidates)
"""

import statistics
from dataclasses import dataclass, field
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, QED, GraphDescriptors

from src.llm.pipeline_layer9 import ValidationResult


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

@dataclass
class MultiFidelityScore:
    smiles:    str
    l1_score:  float   # Layer 9 confidence  [0, 1]
    l2_score:  float   # physicochemical quality  [0, 1]
    l3_score:  float   # synthetic tractability proxy  [0, 1]
    composite: float   # 0.30*L1 + 0.40*L2 + 0.30*L3  [0, 1]
    rank:      int     # 1 = best


@dataclass
class Layer10Result:
    top_candidates: list   # list[MultiFidelityScore], sorted DESC
    screened:       int    # molecules entering L10
    passed_l1:      int    # molecules passing L1 confidence gate
    passed_l2:      int    # molecules passing L2 QED floor
    stats:          dict   # aggregate statistics over all scored molecules


# ---------------------------------------------------------------------------
# Scoring weights and thresholds
# ---------------------------------------------------------------------------

# L1: minimum Layer 9 confidence to enter multi-fidelity pipeline
L1_THRESHOLD  = 0.25

# L2: QED below this floor is too poor for further consideration
L2_QED_FLOOR  = 0.30

# Composite weights (must sum to 1.0)
W_L1 = 0.30
W_L2 = 0.40
W_L3 = 0.30

# LogP sweet-spot [1.0, 3.0] → score 1.0; linear penalty outside
LOGP_IDEAL_LO = 1.0
LOGP_IDEAL_HI = 3.0
LOGP_PENALTY  = 3.0   # normaliser for distance-to-ideal window

# MW sweet-spot centred at 375 Da with ±100 Da half-width
MW_IDEAL      = 375.0
MW_HALF_WIDTH = 100.0

# BertzCT complexity normaliser: 0 = simple, 1000 = very complex
BERTZ_NORM = 1000.0

# Ring count ideal range [2, 4]
RING_IDEAL_LO = 2
RING_IDEAL_HI = 4


# ---------------------------------------------------------------------------
# Per-fidelity scorers (pure functions, each returns a value in [0, 1])
# ---------------------------------------------------------------------------

def _l2_physicochemical(mol) -> float:
    """
    Physicochemical quality score combining:
      - QED (drug-likeness, intrinsically 0–1)
      - LogP proximity to sweet-spot [1–3]
      - BertzCT complexity (lower = simpler synthesis)

    TODO (Phase 2): replace with AutoDock Vina binding energy score
                    after normalising predicted pIC50 to [0, 1].
    """
    qed = QED.qed(mol)
    if qed < L2_QED_FLOOR:
        return 0.0

    logP     = Descriptors.MolLogP(mol)
    logP_gap = max(0.0, max(LOGP_IDEAL_LO - logP, logP - LOGP_IDEAL_HI))
    logP_score = max(0.0, 1.0 - logP_gap / LOGP_PENALTY)

    bertz       = GraphDescriptors.BertzCT(mol)
    bertz_score = max(0.0, 1.0 - bertz / BERTZ_NORM)

    return round((qed + logP_score + bertz_score) / 3.0, 4)


def _l3_tractability(mol) -> float:
    """
    Synthetic tractability proxy combining:
      - MW proximity to 375 Da sweet-spot
      - Ring count reasonableness (2–4 rings is ideal for oral drugs)

    TODO (Phase 3): replace with GROMACS/OpenMM MM-GBSA ΔG_bind
                    after docking pose generation with AutoDock Vina.
    """
    mw       = Descriptors.ExactMolWt(mol)
    mw_score = max(0.0, 1.0 - abs(mw - MW_IDEAL) / MW_HALF_WIDTH)

    n_rings = mol.GetRingInfo().NumRings()
    if RING_IDEAL_LO <= n_rings <= RING_IDEAL_HI:
        ring_score = 1.0
    else:
        ring_dist  = min(
            abs(n_rings - RING_IDEAL_LO),
            abs(n_rings - RING_IDEAL_HI),
        )
        ring_score = max(0.0, 1.0 - ring_dist / 3.0)

    return round((mw_score + ring_score) / 2.0, 4)


# ---------------------------------------------------------------------------
# Main Layer 10 class
# ---------------------------------------------------------------------------

class Layer10MultiFidelityScreening:
    """
    Three-level hierarchical screening of Layer 9 validated candidates.

    Level 1 uses the pre-computed Layer 9 confidence (ML gate).
    Levels 2 and 3 are pure RDKit calculations — no external calls needed.

    The composite score balances immediate drug-likeness (L2) with
    synthetic tractability (L3), weighted by the formal validation result (L1).
    """

    def run(self, validated: list, top_k: int = 10) -> Layer10Result:
        """
        Screen and rank candidates from Layer 9.

        Parameters
        ----------
        validated : All ValidationResult objects from Layer 9 (pass AND fail).
        top_k     : Maximum number of top candidates to return.

        Returns
        -------
        Layer10Result with ranked MultiFidelityScore list.
        """
        screened   = len(validated)
        scored     = []
        passed_l1  = 0
        passed_l2  = 0

        for vr in validated:
            # L1 gate: discard molecules with negligible validation confidence
            if vr.confidence < L1_THRESHOLD:
                continue
            passed_l1 += 1

            mol = Chem.MolFromSmiles(vr.smiles)
            if mol is None:
                continue

            l2 = _l2_physicochemical(mol)
            if l2 <= 0.0:
                continue
            passed_l2 += 1

            l3        = _l3_tractability(mol)
            composite = round(W_L1 * vr.confidence + W_L2 * l2 + W_L3 * l3, 4)

            scored.append(MultiFidelityScore(
                smiles=vr.smiles,
                l1_score=round(vr.confidence, 4),
                l2_score=l2,
                l3_score=l3,
                composite=composite,
                rank=0,
            ))

        # Sort descending by composite and assign ranks
        scored.sort(key=lambda s: s.composite, reverse=True)
        for i, ms in enumerate(scored, start=1):
            ms.rank = i

        top = scored[:top_k]

        composites = [s.composite for s in scored]
        stats: dict = {
            "n_scored":      len(scored),
            "mean_composite": round(statistics.mean(composites), 4) if composites else 0.0,
            "std_composite":  round(statistics.stdev(composites), 4)
                              if len(composites) > 1 else 0.0,
            "max_composite":  round(max(composites), 4) if composites else 0.0,
        }

        return Layer10Result(
            top_candidates=top,
            screened=screened,
            passed_l1=passed_l1,
            passed_l2=passed_l2,
            stats=stats,
        )

    def print_report(self, result: Layer10Result) -> None:
        print("\n[Layer 10] Multi-Fidelity Screening")
        print(f"  {result.screened} in  →  L1 gate: {result.passed_l1}"
              f"  →  L2 gate: {result.passed_l2}"
              f"  →  ranked: {result.stats['n_scored']}")
        print(f"  Composite — mean: {result.stats['mean_composite']:.4f}"
              f"  std: {result.stats['std_composite']:.4f}"
              f"  max: {result.stats['max_composite']:.4f}")
        print()
        print(f"  {'#':>3}  {'L1':>5}  {'L2':>5}  {'L3':>5}  {'Score':>6}  SMILES")
        print("  " + "-" * 78)
        for ms in result.top_candidates:
            abbr = (ms.smiles[:40] + "...") if len(ms.smiles) > 43 else ms.smiles
            print(f"  {ms.rank:>3}  {ms.l1_score:>5.3f}  {ms.l2_score:>5.3f}"
                  f"  {ms.l3_score:>5.3f}  {ms.composite:>6.4f}  {abbr}")


# ---------------------------------------------------------------------------
# Standalone demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from src.llm.pipeline_layer9 import ValidationResult

    demo_results = [
        ValidationResult(
            smiles="CC1=CC=C(C=C1)S(=O)(=O)N", valid=True,
            lipinski_ok=True, veber_ok=True, safety_ok=True,
            proof_verified=True, proof="", confidence=1.0,
            error_type=None, attempts=1,
        ),
        ValidationResult(
            smiles="COc1ccc(CC(N)=O)cc1", valid=True,
            lipinski_ok=True, veber_ok=True, safety_ok=True,
            proof_verified=False, proof="", confidence=0.75,
            error_type="UNRESOLVED_GOALS", attempts=3,
        ),
        ValidationResult(
            smiles="O=C(Nc1cccc(Cl)c1)c1sc2c(c1)CCCC2", valid=True,
            lipinski_ok=True, veber_ok=True, safety_ok=True,
            proof_verified=True, proof="", confidence=1.0,
            error_type=None, attempts=1,
        ),
        ValidationResult(
            smiles="CC(C)(C)C(=O)Nc1sc(CC(N)=O)nc1-c1cccc(F)c1", valid=True,
            lipinski_ok=True, veber_ok=True, safety_ok=False,
            proof_verified=False, proof="", confidence=0.5,
            error_type="NUMERIC_BOUNDS", attempts=2,
        ),
        ValidationResult(
            smiles="COc1ccc(S(=O)(=O)N2CCC(C(N)=O)CC2)cc1", valid=False,
            lipinski_ok=True, veber_ok=True, safety_ok=True,
            proof_verified=False, proof="", confidence=0.75,
            error_type="TYPE_ERROR", attempts=3,
        ),
    ]

    screener = Layer10MultiFidelityScreening()
    result   = screener.run(demo_results, top_k=5)
    screener.print_report(result)
