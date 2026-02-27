"""Layer 6: Target Confidence Scoring.

Aggregates signals from Layers 3–5 into a single normalised confidence
score [0, 1] for each candidate therapeutic target.  High-confidence
targets are forwarded to Layer 7 (molecular generation).

Signal sources and their weights:
┌─────────────────────────────────────────────────────────────────┐
│  Signal                      │ Weight │ Source                  │
├─────────────────────────────────────────────────────────────────┤
│  Open Targets assoc. score   │  0.40  │ Layer 5 (OT GraphQL)    │
│  Betweenness centrality      │  0.25  │ Layer 4 (STRING DB)     │
│  ChEMBL known-drug count     │  0.15  │ Layer 5 (ChEMBL REST)   │
│  Tractability (SM or AB)     │  0.10  │ Layer 5 (OT tractability)│
│  Variant pathogenicity boost │  0.10  │ Layer 3 (VEP / ClinVar) │
├─────────────────────────────────────────────────────────────────┤
│  Penalties                                                       │
│  Essential gene (pLI > 0.9)  │  −0.30 │ Layer 3 (gnomAD)        │
│  OT safety flag present      │  −0.20 │ Layer 5 (OT safety)     │
└─────────────────────────────────────────────────────────────────┘

The weighted sum is clipped to [0, 1] and returned alongside a
human-readable breakdown so the reasoning is fully auditable.

Output feeds into:
  • Layer 7 (molecular generation): chosen target HGNC symbol + EFO disease ID
  • Layer 20 (evidence manager):    per-target confidence scores

Usage:
    from src.pipeline.layer6_target_confidence import Layer6TargetConfidence
    scores = Layer6TargetConfidence().run(l3_result, l4_result, l5_result)
"""

from dataclasses import dataclass, field
from typing import Optional
import math


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

@dataclass
class TargetConfidenceScore:
    gene:               str
    final_score:        float          # [0, 1], weighted aggregate
    ot_component:       float
    centrality_component: float
    chembl_component:   float
    tractability_component: float
    pathogenicity_component: float
    penalty:            float          # negative; from essential/safety flags
    reasoning:          list[str]      = field(default_factory=list)


@dataclass
class Layer6Result:
    scores:         list[TargetConfidenceScore]   # all targets, sorted DESC
    top_target:     Optional[TargetConfidenceScore] = None
    disease_id:     str = ""


# ---------------------------------------------------------------------------
# Scoring weights
# ---------------------------------------------------------------------------

W_OT          = 0.40
W_CENTRALITY  = 0.25
W_CHEMBL      = 0.15
W_TRACT       = 0.10
W_VARIANT     = 0.10

P_ESSENTIAL   = 0.30
P_SAFETY      = 0.20

# ChEMBL drug count normalisation: cap at 20 known drugs → score of 1.0
CHEMBL_CAP = 20


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _log_normalise_chembl(count: int) -> float:
    """Map [0, ∞) drug count to [0, 1] via log scaling capped at CHEMBL_CAP."""
    if count <= 0:
        return 0.0
    return min(math.log1p(count) / math.log1p(CHEMBL_CAP), 1.0)


def _centrality_for_gene(gene: str, l4_result) -> float:
    """Look up betweenness centrality from Layer 4 result."""
    for node in l4_result.nodes:
        if node.gene == gene:
            return min(node.betweenness * 10.0, 1.0)   # scale up; typical values ~0.05–0.15
    return 0.0


def _is_pathogenic_variant(gene: str, l3_result) -> bool:
    """Return True if Layer 3 tagged this gene with a pathogenic consequence."""
    for vc in l3_result.consequences:
        if vc.gene_symbol == gene:
            if vc.impact in ("HIGH", "MODERATE"):
                return True
            if vc.clinical_sig and "pathogenic" in vc.clinical_sig.lower():
                return True
    return False


def _is_essential(gene: str, l3_result) -> bool:
    """Return True if gnomAD marks this gene as loss-of-function intolerant."""
    for gc in l3_result.constraints:
        if gc.gene == gene and not gc.safe_to_target:
            return True
    return False


# ---------------------------------------------------------------------------
# Main Layer 6 class
# ---------------------------------------------------------------------------

class Layer6TargetConfidence:
    """
    Aggregates Layer 3–5 signals into a final ranked confidence list.

    All scoring is deterministic and rule-based (no ML) so the pipeline
    remains fully auditable for regulatory purposes.
    """

    def run(self, l3_result, l4_result, l5_result) -> Layer6Result:
        scores: list[TargetConfidenceScore] = []

        for target in l5_result.ranked_targets:
            gene = target.gene
            reasoning: list[str] = []

            # ── Positive signals ───────────────────────────────────────
            ot_comp = target.ot_association_score * W_OT
            reasoning.append(
                f"OT assoc={target.ot_association_score:.3f} × {W_OT} = {ot_comp:.3f}"
            )

            cent = _centrality_for_gene(gene, l4_result)
            cent_comp = cent * W_CENTRALITY
            reasoning.append(
                f"Centrality={cent:.3f} × {W_CENTRALITY} = {cent_comp:.3f}"
            )

            chembl_norm = _log_normalise_chembl(target.chembl_known_drugs)
            chembl_comp = chembl_norm * W_CHEMBL
            reasoning.append(
                f"ChEMBL drugs={target.chembl_known_drugs} "
                f"(norm={chembl_norm:.3f}) × {W_CHEMBL} = {chembl_comp:.3f}"
            )

            tractable = target.tractable_sm or target.tractable_ab
            tract_comp = float(tractable) * W_TRACT
            reasoning.append(
                f"Tractable={'YES' if tractable else 'NO'} × {W_TRACT} = {tract_comp:.3f}"
            )

            pathogenic = _is_pathogenic_variant(gene, l3_result)
            path_comp  = float(pathogenic) * W_VARIANT
            reasoning.append(
                f"Pathogenic variant={'YES' if pathogenic else 'NO'} "
                f"× {W_VARIANT} = {path_comp:.3f}"
            )

            # ── Penalties ──────────────────────────────────────────────
            penalty = 0.0
            essential = _is_essential(gene, l3_result)
            if essential:
                penalty += P_ESSENTIAL
                reasoning.append(f"PENALTY: essential gene −{P_ESSENTIAL}")

            if target.ot_safety_flags:
                penalty += P_SAFETY
                reasoning.append(
                    f"PENALTY: safety flags {target.ot_safety_flags} −{P_SAFETY}"
                )

            raw_score = ot_comp + cent_comp + chembl_comp + tract_comp + path_comp - penalty
            final = max(0.0, min(1.0, raw_score))
            reasoning.append(f"Final = {final:.4f}")

            scores.append(TargetConfidenceScore(
                gene=gene,
                final_score=final,
                ot_component=ot_comp,
                centrality_component=cent_comp,
                chembl_component=chembl_comp,
                tractability_component=tract_comp,
                pathogenicity_component=path_comp,
                penalty=-penalty,
                reasoning=reasoning,
            ))

        scores.sort(key=lambda s: s.final_score, reverse=True)
        top = scores[0] if scores else None

        return Layer6Result(
            scores=scores,
            top_target=top,
            disease_id=l5_result.chosen_target.ensembl_id if l5_result.chosen_target else "",
        )

    def print_report(self, result: Layer6Result) -> None:
        print("\n[Layer 6] Target Confidence Report")
        print(f"{'Gene':<10} {'Score':>6}  {'OT':>5}  {'Cent':>5}  {'ChEMBL':>6}  Tract  Path  Penalty")
        print("-" * 70)
        for s in result.scores:
            print(
                f"{s.gene:<10} {s.final_score:>6.4f}  "
                f"{s.ot_component:>5.3f}  {s.centrality_component:>5.3f}  "
                f"{s.chembl_component:>6.3f}  {'Y' if s.tractability_component > 0 else 'N':^5}  "
                f"{'Y' if s.pathogenicity_component > 0 else 'N':^4}  "
                f"{s.penalty:>7.3f}"
            )
        if result.top_target:
            print(f"\n→ Recommended target: {result.top_target.gene} "
                  f"(confidence={result.top_target.final_score:.4f})")
            print("  Reasoning:")
            for r in result.top_target.reasoning:
                print(f"    {r}")


# ---------------------------------------------------------------------------
# Quick demo (requires L3–L5 stubs)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Minimal stubs for standalone demo
    from dataclasses import dataclass as dc

    @dc
    class _Node:
        gene: str
        betweenness: float
        degree: int = 0
        pathways: list = field(default_factory=list)

    @dc
    class _L3:
        consequences: list = field(default_factory=list)
        constraints:  list = field(default_factory=list)

    @dc
    class _L4:
        nodes: list = field(default_factory=list)

    from src.pipeline.layer5_target_identification import DrugTarget, Layer5Result

    l3 = _L3()
    l4 = _L4(nodes=[_Node("BCL2", 0.12), _Node("MYC", 0.08)])
    l5 = Layer5Result(ranked_targets=[
        DrugTarget("BCL2", "ENSG00000171791", 0.82, True,  False, 12, True,  False, []),
        DrugTarget("MYC",  "ENSG00000136997", 0.61, False, False,  3, False, False, []),
    ])

    scorer = Layer6TargetConfidence()
    result = scorer.run(l3, l4, l5)
    scorer.print_report(result)
