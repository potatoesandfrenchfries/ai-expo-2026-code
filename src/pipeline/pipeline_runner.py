"""End-to-end drug discovery pipeline orchestrator (Phase 1, Layers 2–10).

Chains all nine implemented layers into a single sequential execution:

    PatientRecord (L2)
      → Layer3GenoPheno       (Ensembl VEP / gnomAD / ClinVar)
      → Layer4CausalPathway   (STRING DB + NetworkX)
      → Layer5TargetID        (Open Targets + ChEMBL)
      → Layer6TargetConfidence (deterministic multi-signal scoring)
      → Layer7MolecularGeneration (Bayesian VAE + ZINC nearest-neighbour)
      → PropertyPredictor (L8) per-candidate
      → Layer9FormalValidator  (Lipinski / Veber / PAINS / Rocq proof)
      → Layer10MultiFidelityScreening (physicochemical + tractability ranking)
      → PipelineResult

Progress, timing, and a unified human-readable report are produced
automatically.  A JSON export is available via DrugDiscoveryPipeline.save_json().

Usage:
    from src.pipeline.pipeline_runner import DrugDiscoveryPipeline, demo_patient
    pipeline = DrugDiscoveryPipeline(n_candidates=30, demo_mode=True)
    result   = pipeline.run(demo_patient())
    pipeline.print_full_report(result)
    pipeline.save_json(result, "output.json")

demo_mode=True:
    Skips real network calls to L3–L5 (Ensembl, STRING, Open Targets, ChEMBL).
    All other layers (L7, L8, L9, L10) run normally.
    Use this for fast local testing without API keys or internet access.
    L9 formal proof generation requires ANTHROPIC_API_KEY; if absent the
    heuristic sub-checks (Lipinski, Veber, PAINS) still run and are reported.
"""

import json
import os
import time
from dataclasses import dataclass, field, asdict
from typing import Optional

from src.pipeline.layer2_patient_intake import PatientRecord
from src.pipeline.layer6_target_confidence import (
    Layer6TargetConfidence, Layer6Result, TargetConfidenceScore,
)
from src.pipeline.layer7_molecular_generation import (
    Layer7MolecularGeneration, Layer7Result,
)
from src.pipeline.layer10_multifidelity_screening import (
    Layer10MultiFidelityScreening, Layer10Result, MultiFidelityScore,
)
from src.llm.pipeline_layer9 import Layer9FormalValidator, ValidationResult
from src.llm.feedback_controller import FeedbackController


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class PipelineResult:
    patient_id:         str
    disease:            str
    top_target:         str
    target_confidence:  float
    top_candidates:     list   # list[MultiFidelityScore]
    validation_summary: dict   # {total, passed, avg_confidence, proof_rate}
    layer_timing_s:     dict   # {L3: 4.2, L4: 1.1, ...}
    success:            bool
    error_layer:        Optional[str] = None
    error_message:      Optional[str] = None


# ---------------------------------------------------------------------------
# Demo / stub patient
# ---------------------------------------------------------------------------

def demo_patient() -> PatientRecord:
    """
    Return the canonical PT-0001 DLBCL demonstration patient.
    This is the same example used in layer2_patient_intake.py so
    all downstream layers produce deterministic stub data.
    """
    from src.pipeline.layer2_patient_intake import GeneVariant, LabValue
    return PatientRecord(
        patient_id="PT-0001",
        age=54,
        sex="F",
        presenting_symptoms=["fatigue", "unexplained weight loss", "night sweats"],
        confirmed_diagnosis="C83.3 Diffuse large B-cell lymphoma",
        disease_subtype="GCB",
        severity="severe",
        genes_of_interest=["BCL2", "MYC", "TP53"],
        gene_variants=[
            GeneVariant(
                gene="BCL2",
                variant_id="NM_000633.3:c.118A>G",
                zygosity="heterozygous",
                clinical_significance="pathogenic",
            ),
        ],
        lab_values=[
            LabValue(
                test_name="LDH",
                value=580.0,
                unit="U/L",
                reference_range="135–225 U/L",
                flag="high",
            ),
        ],
        notes="Referred from oncology. Awaiting PET-CT.",
    )


# ---------------------------------------------------------------------------
# Demo stubs for L3–L5 (no network calls)
# ---------------------------------------------------------------------------

def _make_demo_l3(patient: PatientRecord):
    """Minimal Layer3Result compatible with Layer6TargetConfidence."""
    from src.pipeline.layer3_genotype_phenotype import (
        Layer3Result, VariantConsequence, GeneConstraint,
    )
    consequences = []
    for gv in patient.gene_variants:
        consequences.append(VariantConsequence(
            variant_id=gv.variant_id,
            gene_symbol=gv.gene,
            transcript_id="ENST_STUB",
            consequence="missense_variant",
            impact="HIGH" if gv.clinical_significance == "pathogenic" else "MODERATE",
            sift="deleterious",
            polyphen="probably_damaging",
            clinical_sig=gv.clinical_significance,
        ))
    constraints = [
        GeneConstraint(gene="BCL2", pLI=0.02, LOEUF=0.82, safe_to_target=True),
        GeneConstraint(gene="MYC",  pLI=0.10, LOEUF=0.91, safe_to_target=True),
        GeneConstraint(gene="TP53", pLI=0.01, LOEUF=1.10, safe_to_target=True),
    ]
    dysfunctional = list({gv.gene for gv in patient.gene_variants})
    return Layer3Result(
        consequences=consequences,
        constraints=constraints,
        dysfunctional_proteins=dysfunctional,
    )


def _make_demo_l4(dysfunctional_proteins: list):
    """Minimal Layer4Result compatible with Layer6TargetConfidence."""
    from src.pipeline.layer4_causal_pathway import Layer4Result, ProteinNode
    import networkx as nx

    nodes = [
        ProteinNode(gene="BCL2", degree=14, betweenness=0.12,
                    pathways=["Apoptosis", "BCL-2 family proteins"]),
        ProteinNode(gene="MYC",  degree=9,  betweenness=0.08,
                    pathways=["MYC targets", "Cell cycle"]),
        ProteinNode(gene="TP53", degree=18, betweenness=0.15,
                    pathways=["TP53 regulation", "Apoptosis"]),
    ]
    # Add any dysfunctional proteins not already in stub nodes
    stub_genes = {n.gene for n in nodes}
    for gene in dysfunctional_proteins:
        if gene not in stub_genes:
            nodes.append(ProteinNode(gene=gene, degree=5, betweenness=0.05))

    candidate_targets = [n.gene for n in sorted(nodes, key=lambda n: -n.betweenness)]

    G = nx.DiGraph()
    for n in nodes:
        G.add_node(n.gene)
    for i, a in enumerate(nodes):
        for b in nodes[i + 1:]:
            G.add_edge(a.gene, b.gene, weight=0.7)

    return Layer4Result(
        graph=G,
        nodes=nodes,
        candidate_targets=candidate_targets,
        causal_chains={p: [p, "BCL2"] for p in dysfunctional_proteins},
    )


def _make_demo_l5(candidate_targets: list):
    """Minimal Layer5Result compatible with Layer6TargetConfidence."""
    from src.pipeline.layer5_target_identification import Layer5Result, DrugTarget

    STUB_DATA = {
        "BCL2": DrugTarget(
            gene="BCL2", ensembl_id="ENSG00000171791",
            ot_association_score=0.82, tractable_sm=True, tractable_ab=False,
            chembl_known_drugs=12, has_binding_pocket=True,
            essential_gene=False, ot_safety_flags=[],
        ),
        "MYC": DrugTarget(
            gene="MYC", ensembl_id="ENSG00000136997",
            ot_association_score=0.61, tractable_sm=False, tractable_ab=False,
            chembl_known_drugs=3, has_binding_pocket=False,
            essential_gene=False, ot_safety_flags=[],
        ),
        "TP53": DrugTarget(
            gene="TP53", ensembl_id="ENSG00000141510",
            ot_association_score=0.71, tractable_sm=True, tractable_ab=False,
            chembl_known_drugs=6, has_binding_pocket=True,
            essential_gene=False, ot_safety_flags=[],
        ),
    }

    ranked = []
    for gene in candidate_targets:
        if gene in STUB_DATA:
            ranked.append(STUB_DATA[gene])
        else:
            # Generic fallback for unknown genes
            ranked.append(DrugTarget(
                gene=gene, ensembl_id="ENSG_STUB",
                ot_association_score=0.30, tractable_sm=False, tractable_ab=False,
                chembl_known_drugs=0, has_binding_pocket=False,
                essential_gene=False, ot_safety_flags=[],
            ))

    ranked.sort(key=lambda t: t.ot_association_score, reverse=True)
    chosen = next((t for t in ranked if t.tractable_sm or t.tractable_ab), ranked[0] if ranked else None)

    return Layer5Result(ranked_targets=ranked, chosen_target=chosen)


# ---------------------------------------------------------------------------
# L8 property computation (RDKit direct — avoids training dependency)
# ---------------------------------------------------------------------------

def _compute_l8_properties(smiles_list: list) -> dict:
    """
    Compute L8 properties (QED, LogP, MW, TPSA, HBD, HBA) directly with
    RDKit for each candidate SMILES.  Returns a dict keyed by SMILES.

    This mirrors the output format of property_predictor.predict_smiles()
    but uses deterministic RDKit values instead of the trained MLP.
    Avoids the data/property_predictor.pt checkpoint dependency at runtime.
    """
    from rdkit.Chem import Descriptors, rdMolDescriptors, QED

    results = {}
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        results[smi] = {
            "QED":  {"mean": round(QED.qed(mol), 4),                            "std": 0.0},
            "LogP": {"mean": round(Descriptors.MolLogP(mol), 4),                "std": 0.0},
            "MW":   {"mean": round(Descriptors.ExactMolWt(mol), 4),             "std": 0.0},
            "TPSA": {"mean": round(Descriptors.TPSA(mol), 4),                   "std": 0.0},
            "HBD":  {"mean": float(rdMolDescriptors.CalcNumHBD(mol)),           "std": 0.0},
            "HBA":  {"mean": float(rdMolDescriptors.CalcNumHBA(mol)),           "std": 0.0},
        }
    return results


# ---------------------------------------------------------------------------
# Main pipeline class
# ---------------------------------------------------------------------------

class DrugDiscoveryPipeline:
    """
    End-to-end Phase 1 drug discovery pipeline.

    Parameters
    ----------
    n_candidates : Number of molecular candidates to generate in Layer 7.
    demo_mode    : Skip network calls to L3–L5 and use stub data instead.
                   Useful for rapid local testing without API access.
    """

    def __init__(self, n_candidates: int = 50, demo_mode: bool = False):
        self.n_candidates = n_candidates
        self.demo_mode    = demo_mode

    # -----------------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------------

    def run(self, patient_record: PatientRecord) -> PipelineResult:
        """Execute the full L2–L10 pipeline for a given patient record."""
        self._print_banner(patient_record)

        timing: dict = {}
        l3_result = l4_result = l5_result = l6_result = None
        l7_result = l8_props  = l9_results = l10_result = None

        # ── Layer 3: Genotype → Phenotype ───────────────────────────────────
        t0 = time.time()
        try:
            if self.demo_mode:
                l3_result = _make_demo_l3(patient_record)
            else:
                from src.pipeline.layer3_genotype_phenotype import Layer3GenoPheno
                l3_result = Layer3GenoPheno().run(
                    patient_record.hgvs_notations(),
                    patient_record.genes_of_interest,
                )
            timing["L3"] = round(time.time() - t0, 2)
            self._step_ok("L3", "Genotype-phenotype",
                          f"{len(l3_result.consequences)} consequences"
                          f" | {len(l3_result.dysfunctional_proteins)} dysfunctional genes")
        except Exception as exc:
            return self._abort("L3", exc, timing, patient_record)

        # ── Layer 4: Causal Pathway ──────────────────────────────────────────
        t0 = time.time()
        try:
            if self.demo_mode:
                l4_result = _make_demo_l4(l3_result.dysfunctional_proteins)
            else:
                from src.pipeline.layer4_causal_pathway import Layer4CausalPathway
                safe_genes = [gc.gene for gc in l3_result.constraints if gc.safe_to_target]
                l4_result  = Layer4CausalPathway().run(
                    l3_result.dysfunctional_proteins,
                    safe_genes=safe_genes,
                )
            timing["L4"] = round(time.time() - t0, 2)
            n_edges = l4_result.graph.number_of_edges() if hasattr(l4_result.graph, "number_of_edges") else "?"
            self._step_ok("L4", "Causal pathways",
                          f"{len(l4_result.nodes)} proteins | {n_edges} interactions")
        except Exception as exc:
            return self._abort("L4", exc, timing, patient_record)

        # ── Layer 5: Target Identification ───────────────────────────────────
        t0 = time.time()
        try:
            if self.demo_mode:
                l5_result = _make_demo_l5(l4_result.candidate_targets)
            else:
                from src.pipeline.layer5_target_identification import Layer5TargetID
                l5_result = Layer5TargetID().run(l4_result.candidate_targets)
            timing["L5"] = round(time.time() - t0, 2)
            self._step_ok("L5", "Target identification",
                          f"{len(l5_result.ranked_targets)} candidates")
        except Exception as exc:
            return self._abort("L5", exc, timing, patient_record)

        # ── Layer 6: Target Confidence Scoring ──────────────────────────────
        t0 = time.time()
        try:
            scorer    = Layer6TargetConfidence()
            l6_result = scorer.run(l3_result, l4_result, l5_result)
            timing["L6"] = round(time.time() - t0, 2)
            top = l6_result.top_target
            self._step_ok("L6", "Target confidence",
                          f"{top.gene}  score={top.final_score:.4f}" if top else "no target")
        except Exception as exc:
            return self._abort("L6", exc, timing, patient_record)

        # ── Layer 7: Molecular Generation ───────────────────────────────────
        t0 = time.time()
        feedback = FeedbackController()
        try:
            layer7    = Layer7MolecularGeneration()
            l7_result = layer7.run(l6_result, n_candidates=self.n_candidates)
            timing["L7"] = round(time.time() - t0, 2)
            self._step_ok("L7", "Molecular generation",
                          f"{l7_result.n_valid} candidates generated")
        except Exception as exc:
            return self._abort("L7", exc, timing, patient_record)

        # ── Layer 8: Property Prediction ─────────────────────────────────────
        t0 = time.time()
        try:
            smiles_list = [c.smiles for c in l7_result.candidates]
            l8_props    = _compute_l8_properties(smiles_list)
            timing["L8"] = round(time.time() - t0, 2)
            self._step_ok("L8", "Property prediction",
                          f"{len(l8_props)}/{len(smiles_list)} computed")
        except Exception as exc:
            return self._abort("L8", exc, timing, patient_record)

        # ── Layer 9: Formal Validation ────────────────────────────────────────
        t0 = time.time()
        try:
            validator  = Layer9FormalValidator()
            l9_results = validator.validate_batch(smiles_list)
            timing["L9"] = round(time.time() - t0, 2)

            n_pass = sum(1 for r in l9_results if r.valid)
            n_fail = len(l9_results) - n_pass
            self._step_ok("L9", "Formal validation",
                          f"{n_pass} pass | {n_fail} fail")

            # Feed failures back into L7 constraints for future refinement
            failures = [r for r in l9_results if not r.valid]
            if failures:
                feedback.extract_constraints_from_failures(failures)
        except Exception as exc:
            return self._abort("L9", exc, timing, patient_record)

        # ── Layer 10: Multi-Fidelity Screening ───────────────────────────────
        t0 = time.time()
        try:
            screener   = Layer10MultiFidelityScreening()
            l10_result = screener.run(l9_results, top_k=min(10, self.n_candidates))
            timing["L10"] = round(time.time() - t0, 2)
            self._step_ok("L10", "Multi-fidelity screen",
                          f"{len(l10_result.top_candidates)} top candidates selected")
        except Exception as exc:
            return self._abort("L10", exc, timing, patient_record)

        # ── Aggregate result ─────────────────────────────────────────────────
        avg_conf = (
            round(sum(r.confidence for r in l9_results) / len(l9_results), 4)
            if l9_results else 0.0
        )
        proof_rate = (
            round(sum(1 for r in l9_results if r.proof_verified) / len(l9_results), 4)
            if l9_results else 0.0
        )

        return PipelineResult(
            patient_id=patient_record.patient_id,
            disease=patient_record.confirmed_diagnosis,
            top_target=l6_result.top_target.gene if l6_result.top_target else "UNKNOWN",
            target_confidence=l6_result.top_target.final_score if l6_result.top_target else 0.0,
            top_candidates=l10_result.top_candidates,
            validation_summary={
                "total":      len(l9_results),
                "passed":     sum(1 for r in l9_results if r.valid),
                "avg_confidence": avg_conf,
                "proof_rate": proof_rate,
            },
            layer_timing_s=timing,
            success=True,
        )

    def print_full_report(self, result: PipelineResult) -> None:
        """Print the unified pipeline summary report."""
        total_s = round(sum(result.layer_timing_s.values()), 1)

        print()
        print("=" * 62)
        print(f"  TARGET        : {result.top_target}")
        print(f"  CONFIDENCE    : {result.target_confidence:.4f}")
        print(f"  PATIENT       : {result.patient_id}  ({result.disease})")
        print("=" * 62)

        if not result.success:
            print(f"\n  PIPELINE FAILED at {result.error_layer}: {result.error_message}")
            return

        val = result.validation_summary
        print(f"\n  Validation    : {val.get('passed', 0)}/{val.get('total', 0)} passed"
              f"  avg conf={val.get('avg_confidence', 0):.3f}"
              f"  proof rate={val.get('proof_rate', 0):.2f}")

        print(f"\n  {'#':>3}  {'Score':>6}  {'QED?':>5}  {'MW':>6}  {'LogP':>5}  SMILES")
        print("  " + "-" * 68)

        from rdkit.Chem import Descriptors, QED
        for ms in result.top_candidates:
            mol = Chem.MolFromSmiles(ms.smiles)
            qed_v = round(QED.qed(mol), 2)     if mol else "N/A"
            mw_v  = round(Descriptors.ExactMolWt(mol), 1) if mol else "N/A"
            lp_v  = round(Descriptors.MolLogP(mol), 2)    if mol else "N/A"
            abbr  = (ms.smiles[:35] + "...") if len(ms.smiles) > 38 else ms.smiles
            print(f"  {ms.rank:>3}  {ms.composite:>6.4f}  {str(qed_v):>5}  "
                  f"{str(mw_v):>6}  {str(lp_v):>5}  {abbr}")

        print()
        print(f"  Completed in {total_s}s")
        if result.error_layer:
            print(f"  WARNING: pipeline recovered from error at {result.error_layer}:"
                  f" {result.error_message}")

    def save_json(self, result: PipelineResult, path: str) -> None:
        """Serialise PipelineResult to a JSON file."""
        data = {
            "patient_id":         result.patient_id,
            "disease":            result.disease,
            "target":             {
                "gene":       result.top_target,
                "confidence": result.target_confidence,
            },
            "candidates": [
                {
                    "rank":      ms.rank,
                    "smiles":    ms.smiles,
                    "composite": ms.composite,
                    "l1_score":  ms.l1_score,
                    "l2_score":  ms.l2_score,
                    "l3_score":  ms.l3_score,
                }
                for ms in result.top_candidates
            ],
            "validation_summary": result.validation_summary,
            "layer_timing_s":     result.layer_timing_s,
            "success":            result.success,
            "error_layer":        result.error_layer,
            "error_message":      result.error_message,
        }
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
        print(f"\n  Results saved → {path}")

    # -----------------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------------

    @staticmethod
    def _print_banner(patient: PatientRecord) -> None:
        print()
        print("╔" + "═" * 60 + "╗")
        print("║  KNOTWORKING DRUG DISCOVERY PIPELINE  —  Phase 1" + " " * 10 + "║")
        print("╚" + "═" * 60 + "╝")
        print(f"  Patient  : {patient.patient_id}  |  Age {patient.age}  |  Sex {patient.sex}")
        print(f"  Disease  : {patient.confirmed_diagnosis}")
        print(f"  Genes    : {', '.join(patient.genes_of_interest)}")
        print()

    @staticmethod
    def _step_ok(tag: str, label: str, detail: str) -> None:
        col_width = 26
        padded = f"[{tag}] {label}"
        dots   = "." * max(1, col_width - len(padded))
        print(f"  {padded}{dots}  ✓  {detail}")

    @staticmethod
    def _abort(layer: str, exc: Exception, timing: dict, patient: PatientRecord) -> PipelineResult:
        print(f"\n  [{layer}] ✗  ERROR: {exc}")
        return PipelineResult(
            patient_id=patient.patient_id,
            disease=patient.confirmed_diagnosis,
            top_target="UNKNOWN",
            target_confidence=0.0,
            top_candidates=[],
            validation_summary={},
            layer_timing_s=timing,
            success=False,
            error_layer=layer,
            error_message=str(exc),
        )
