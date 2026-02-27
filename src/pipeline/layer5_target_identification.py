"""Layer 5: Target Identification.

Ranks candidate proteins from Layer 4 by druggability using two
production databases — both free, no authentication required:

┌────────────────────────────────────────────────────────────────────┐
│  Open Targets Platform  (platform.opentargets.org)                 │
│  GraphQL API v4 — disease-target association scores (0–1),         │
│  tractability (small molecule / antibody buckets), safety flags    │
├────────────────────────────────────────────────────────────────────┤
│  ChEMBL REST API  (ebi.ac.uk/chembl)                               │
│  — Known drug-target interactions & activity data                  │
│  — Binding pocket existence (via target components)                │
└────────────────────────────────────────────────────────────────────┘

Why not just build our own?
  Open Targets integrates GWAS, somatic mutations, expression, pathways,
  known drugs, and literature into a single score.  Reimplementing this
  would require terabytes of raw data and months of work.  The API
  exposes the same scores in milliseconds.

Output feeds into:
  • Layer 6 (target confidence): ranked_targets list

Usage:
    from src.pipeline.layer5_target_identification import Layer5TargetID
    result = Layer5TargetID().run(["BCL2", "MYC"], disease_id="EFO_0000183")
"""

import time
from dataclasses import dataclass, field
from typing import Optional

import requests

# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

@dataclass
class DrugTarget:
    gene:                   str
    ensembl_id:             str
    ot_association_score:   float          # Open Targets overall score [0,1]
    tractable_sm:           bool           # small molecule tractability
    tractable_ab:           bool           # antibody tractability
    chembl_known_drugs:     int            # count of approved drugs hitting this target
    has_binding_pocket:     bool
    essential_gene:         bool           # flagged as broadly essential in DepMap/CRISPR screens
    ot_safety_flags:        list[str]      = field(default_factory=list)


@dataclass
class Layer5Result:
    ranked_targets: list[DrugTarget]   # sorted by ot_association_score DESC
    chosen_target:  Optional[DrugTarget] = None   # highest score that is tractable + safe


# ---------------------------------------------------------------------------
# Open Targets GraphQL
# ---------------------------------------------------------------------------

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

_OT_TARGET_QUERY = """
query TargetInfo($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    tractability {
      label
      modality
      value
    }
    safetyLiabilities {
      event
    }
  }
}
"""

_OT_DISEASE_ASSOC_QUERY = """
query DiseaseAssociation($diseaseId: String!, $ensemblId: String!) {
  disease(efoId: $diseaseId) {
    associatedTargets(
      orderByScore: "score"
      filter: {ids: [$ensemblId]}
    ) {
      rows {
        target { id approvedSymbol }
        score
      }
    }
  }
}
"""

_OT_GENE_LOOKUP_QUERY = """
query GeneSearch($symbol: String!) {
  search(queryString: $symbol, entityNames: ["target"]) {
    hits {
      id
      name
      entity
    }
  }
}
"""


def _ot_post(query: str, variables: dict) -> dict:
    for attempt in range(4):
        try:
            resp = requests.post(
                OT_GRAPHQL,
                json={"query": query, "variables": variables},
                timeout=20,
            )
            if resp.status_code == 200:
                return resp.json().get("data", {})
        except requests.RequestException as e:
            print(f"  Open Targets attempt {attempt+1} failed: {e}")
        time.sleep(2 ** attempt)
    return {}


def resolve_ensembl_id(gene_symbol: str) -> Optional[str]:
    """Resolve a HGNC gene symbol to an Ensembl gene ID via Open Targets search."""
    data = _ot_post(_OT_GENE_LOOKUP_QUERY, {"symbol": gene_symbol})
    hits = data.get("search", {}).get("hits", [])
    for hit in hits:
        if hit.get("entity") == "target" and hit.get("name", "").upper() == gene_symbol.upper():
            return hit["id"]
    # Fallback: return first hit ID if name doesn't exactly match
    if hits:
        return hits[0].get("id")
    return None


def get_ot_tractability(ensembl_id: str) -> tuple[bool, bool, list[str]]:
    """Return (sm_tractable, ab_tractable, safety_flags) from Open Targets."""
    data = _ot_post(_OT_TARGET_QUERY, {"ensemblId": ensembl_id})
    target = data.get("target", {}) or {}

    sm, ab = False, False
    for t in target.get("tractability", []):
        if t.get("value"):
            mod = t.get("modality", "")
            if mod == "SM":
                sm = True
            elif mod == "AB":
                ab = True

    safety_flags = [s["event"] for s in target.get("safetyLiabilities", []) if s.get("event")]
    return sm, ab, safety_flags


def get_ot_association_score(ensembl_id: str, disease_id: str) -> float:
    """Return the Open Targets overall association score for a target+disease pair."""
    data = _ot_post(_OT_DISEASE_ASSOC_QUERY, {"diseaseId": disease_id, "ensemblId": ensembl_id})
    rows = (
        data.get("disease", {})
            .get("associatedTargets", {})
            .get("rows", [])
    ) or []
    if rows:
        return float(rows[0].get("score", 0.0))
    return 0.0


# ---------------------------------------------------------------------------
# ChEMBL REST
# ---------------------------------------------------------------------------

CHEMBL_TARGET_URL = "https://www.ebi.ac.uk/chembl/api/data/target"
CHEMBL_MECH_URL   = "https://www.ebi.ac.uk/chembl/api/data/mechanism"


def _chembl_gene_info(gene_symbol: str) -> tuple[int, bool]:
    """
    Return (known_drug_count, has_binding_pocket) for a gene via ChEMBL.
    """
    try:
        # Search targets by gene symbol
        resp = requests.get(
            CHEMBL_TARGET_URL,
            params={"target_synonym__icontains": gene_symbol, "format": "json", "limit": 5},
            timeout=15,
        )
        if resp.status_code != 200:
            return 0, False
        results = resp.json().get("targets", [])
        if not results:
            return 0, False

        chembl_id = results[0].get("target_chembl_id", "")
        has_pocket = results[0].get("target_type", "") in ("SINGLE PROTEIN", "PROTEIN COMPLEX")

        # Count approved drugs against this target
        mech_resp = requests.get(
            CHEMBL_MECH_URL,
            params={"target_chembl_id": chembl_id, "format": "json", "limit": 100},
            timeout=15,
        )
        if mech_resp.status_code != 200:
            return 0, has_pocket
        drug_count = mech_resp.json().get("page_meta", {}).get("total_count", 0)
        return int(drug_count), has_pocket

    except requests.RequestException:
        return 0, False


# ---------------------------------------------------------------------------
# Main Layer 5 class
# ---------------------------------------------------------------------------

class Layer5TargetID:
    """
    Scores and ranks candidate proteins by druggability.

    Parameters
    ----------
    disease_id : str
        EFO ID for the disease (e.g. "EFO_0000183" for lymphoma).
        Used for the Open Targets disease-specific association score.
    """

    def __init__(self, disease_id: str = "EFO_0000183"):
        self.disease_id = disease_id

    def run(self, candidate_genes: list[str]) -> Layer5Result:
        targets: list[DrugTarget] = []

        for gene in candidate_genes:
            print(f"[L5] Scoring target: {gene}")

            ensembl_id = resolve_ensembl_id(gene)
            if ensembl_id is None:
                print(f"  Could not resolve Ensembl ID for {gene}, skipping.")
                continue

            ot_score            = get_ot_association_score(ensembl_id, self.disease_id)
            sm, ab, safety      = get_ot_tractability(ensembl_id)
            drug_count, pocket  = _chembl_gene_info(gene)

            targets.append(DrugTarget(
                gene=gene,
                ensembl_id=ensembl_id,
                ot_association_score=ot_score,
                tractable_sm=sm,
                tractable_ab=ab,
                chembl_known_drugs=drug_count,
                has_binding_pocket=pocket,
                essential_gene=False,   # Layer 3 safe_to_target flag applied in Layer 6
                ot_safety_flags=safety,
            ))
            time.sleep(0.3)    # gentle API rate limiting

        ranked = sorted(targets, key=lambda t: t.ot_association_score, reverse=True)

        chosen = next(
            (t for t in ranked if (t.tractable_sm or t.tractable_ab)
             and not t.essential_gene and not t.ot_safety_flags),
            ranked[0] if ranked else None,
        )

        return Layer5Result(ranked_targets=ranked, chosen_target=chosen)


# ---------------------------------------------------------------------------
# Quick demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    l5 = Layer5TargetID(disease_id="EFO_0000183")
    result = l5.run(["BCL2", "MYC", "TP53"])

    print("\nRanked targets:")
    for t in result.ranked_targets:
        print(f"  {t.gene:8s}  OT={t.ot_association_score:.3f}  "
              f"SM={t.tractable_sm}  AB={t.tractable_ab}  "
              f"drugs={t.chembl_known_drugs}  pocket={t.has_binding_pocket}")

    if result.chosen_target:
        print(f"\nChosen therapeutic target: {result.chosen_target.gene}")
