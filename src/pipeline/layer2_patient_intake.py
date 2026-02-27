"""Layer 2: Patient data intake schema.

L2 is NOT a machine-learning layer.  A clinician enters structured patient
data (symptoms, diagnosis, genes of interest, lab values) which becomes the
seed for all downstream computational layers.  The schema is intentionally
minimal — just enough to identify a therapeutic target.

Downstream consumers:
  • Layer 3 (genotype-phenotype): reads gene_variants
  • Layer 4 (causal pathway):     reads genes_of_interest
  • Layer 5 (target ID):          reads confirmed_diagnosis + genes_of_interest

Usage:
    from src.pipeline.layer2_patient_intake import PatientRecord, load_record, save_record
"""

import json
import os
from dataclasses import dataclass, field, asdict
from typing import Optional


# ---------------------------------------------------------------------------
# Sub-schemas
# ---------------------------------------------------------------------------

@dataclass
class GeneVariant:
    """A single genetic variant reported for the patient."""
    gene: str                    # HGNC symbol, e.g. "BRCA1"
    variant_id: str              # dbSNP rs ID or HGVS notation
    zygosity: str                # "heterozygous" | "homozygous" | "unknown"
    clinical_significance: str   # "pathogenic" | "likely_pathogenic" | "VUS" | "benign"


@dataclass
class LabValue:
    """A single quantitative lab measurement."""
    test_name: str
    value: float
    unit: str
    reference_range: str         # e.g. "3.5–5.0 mmol/L"
    flag: str                    # "normal" | "high" | "low" | "critical"


# ---------------------------------------------------------------------------
# Core patient record
# ---------------------------------------------------------------------------

@dataclass
class PatientRecord:
    """
    Structured patient data provided by the attending clinician.

    Fields the pipeline uses directly:
        confirmed_diagnosis   → search string for L5 Open Targets query
        genes_of_interest     → protein list for L4 STRING DB query
        gene_variants         → HGVS notations for L3 Ensembl VEP query
    """

    # --- Identity (anonymised) ---
    patient_id:            str
    age:                   int
    sex:                   str                       # "M" | "F" | "other"

    # --- Clinical picture ---
    presenting_symptoms:   list[str]                 # free-text symptom list
    confirmed_diagnosis:   str                       # ICD-10 code + plain text, e.g. "C50 Breast cancer"
    disease_subtype:       Optional[str] = None      # e.g. "HER2-positive"
    severity:              Optional[str] = None      # "mild" | "moderate" | "severe"

    # --- Genomics (may be empty if WGS not available) ---
    genes_of_interest:     list[str]     = field(default_factory=list)   # HGNC symbols
    gene_variants:         list[GeneVariant] = field(default_factory=list)

    # --- Lab results ---
    lab_values:            list[LabValue]    = field(default_factory=list)

    # --- Clinician notes ---
    notes:                 Optional[str] = None

    # ---------------------------------------------------------------------------
    # Validation
    # ---------------------------------------------------------------------------

    def validate(self) -> list[str]:
        """Return a list of validation warnings (empty = clean)."""
        warnings: list[str] = []
        if not self.patient_id:
            warnings.append("patient_id is empty.")
        if not self.confirmed_diagnosis:
            warnings.append("confirmed_diagnosis is required for downstream layers.")
        if self.age <= 0 or self.age > 130:
            warnings.append(f"Suspicious age: {self.age}")
        for gv in self.gene_variants:
            if not gv.variant_id:
                warnings.append(f"Gene variant for {gv.gene} has no variant_id.")
        return warnings

    # ---------------------------------------------------------------------------
    # L3 helper: extract HGVS notations for Ensembl VEP
    # ---------------------------------------------------------------------------

    def hgvs_notations(self) -> list[str]:
        """Return all HGVS notations ready for Ensembl VEP POST."""
        return [gv.variant_id for gv in self.gene_variants if gv.variant_id]

    # ---------------------------------------------------------------------------
    # Serialisation
    # ---------------------------------------------------------------------------

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "PatientRecord":
        variants = [GeneVariant(**v) for v in data.pop("gene_variants", [])]
        labs     = [LabValue(**l)    for l in data.pop("lab_values", [])]
        return cls(**data, gene_variants=variants, lab_values=labs)


# ---------------------------------------------------------------------------
# Persistence helpers
# ---------------------------------------------------------------------------

def save_record(record: PatientRecord, path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump(record.to_dict(), f, indent=2)
    print(f"Patient record saved → {path}")


def load_record(path: str) -> PatientRecord:
    with open(path) as f:
        return PatientRecord.from_dict(json.load(f))


# ---------------------------------------------------------------------------
# Demo / quick test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    demo = PatientRecord(
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
            )
        ],
        lab_values=[
            LabValue(
                test_name="LDH",
                value=580.0,
                unit="U/L",
                reference_range="135–225 U/L",
                flag="high",
            )
        ],
        notes="Referred from oncology. Awaiting PET-CT.",
    )

    warnings = demo.validate()
    if warnings:
        print("Validation warnings:", warnings)
    else:
        print("Record valid.")

    save_record(demo, "data/demo_patient.json")
    loaded = load_record("data/demo_patient.json")
    print("HGVS notations for L3:", loaded.hgvs_notations())
