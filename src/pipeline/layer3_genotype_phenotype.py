"""Layer 3: Genotype-to-Phenotype Association.

Wraps three existing, production-grade bioinformatics services:

┌──────────────────────────────────────────────────────────────────┐
│  Tool               │ What it does                               │
├──────────────────────────────────────────────────────────────────┤
│  Ensembl VEP REST   │ Variant → transcript consequence + impact  │
│  (EMBL-EBI)         │ No auth required. Rate: 200 variants/POST. │
├──────────────────────────────────────────────────────────────────┤
│  gnomAD GraphQL     │ Population allele frequency + constraint   │
│  (Broad Institute)  │ No auth required.                          │
├──────────────────────────────────────────────────────────────────┤
│  ClinVar E-utilities│ Clinical significance for known variants   │
│  (NCBI)             │ No auth required (email recommended).      │
└──────────────────────────────────────────────────────────────────┘

Why these three and not something custom?
- Ensembl VEP is the gold standard for variant annotation; it underlies
  AlphaMissense, SIFT, PolyPhen-2, CADD, and most clinical pipelines.
- gnomAD provides constraint scores (pLI, LOEUF) that tell us whether
  a gene tolerates loss-of-function — critical for L5 target safety.
- ClinVar gives direct clinical labels, avoiding the need for ML here.

No ML is used in this layer. All annotation is database/API-driven.

Output feeds into:
  • Layer 4 (causal pathway): dysfunctional_proteins
  • Layer 5 (target ID):      constrained_genes (pLI > 0.9 → avoid)

Usage:
    from src.pipeline.layer3_genotype_phenotype import Layer3GenoPheno
    result = Layer3GenoPheno().run(["NM_000633.3:c.118A>G"], genes=["BCL2"])
"""

import time
from dataclasses import dataclass, field
from typing import Optional

import requests

# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

@dataclass
class VariantConsequence:
    variant_id:    str
    gene_symbol:   str
    transcript_id: str
    consequence:   str     # e.g. "missense_variant"
    impact:        str     # "HIGH" | "MODERATE" | "LOW" | "MODIFIER"
    sift:          Optional[str] = None      # "deleterious" | "tolerated"
    polyphen:      Optional[str] = None      # "probably_damaging" etc.
    clinical_sig:  Optional[str] = None      # from ClinVar


@dataclass
class GeneConstraint:
    gene:  str
    pLI:   float    # prob. loss-of-function intolerant (gnomAD); >0.9 → avoid targeting
    LOEUF: float    # loss-of-function observed/expected upper bound; <0.35 → constrained
    safe_to_target: bool  # True if gene can be disrupted without likely lethality


@dataclass
class Layer3Result:
    consequences:  list[VariantConsequence] = field(default_factory=list)
    constraints:   list[GeneConstraint]     = field(default_factory=list)
    dysfunctional_proteins: list[str]       = field(default_factory=list)


# ---------------------------------------------------------------------------
# Ensembl VEP
# ---------------------------------------------------------------------------

VEP_URL = "https://rest.ensembl.org/vep/human/hgvs"
_VEP_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


def _query_vep(hgvs_list: list[str]) -> list[dict]:
    """POST up to 200 HGVS notations to Ensembl VEP and return raw JSON."""
    if not hgvs_list:
        return []

    payload = {"hgvs_notations": hgvs_list}
    for attempt in range(4):
        try:
            resp = requests.post(VEP_URL, json=payload, headers=_VEP_HEADERS, timeout=30)
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 429:
                wait = int(resp.headers.get("Retry-After", 5))
                time.sleep(wait)
                continue
        except requests.RequestException as e:
            print(f"  VEP attempt {attempt+1} failed: {e}")
        time.sleep(2 ** attempt)
    return []


def parse_vep(raw: list[dict]) -> list[VariantConsequence]:
    results: list[VariantConsequence] = []
    for entry in raw:
        variant_id = entry.get("id", entry.get("input", "?"))
        for tc in entry.get("transcript_consequences", []):
            results.append(VariantConsequence(
                variant_id=variant_id,
                gene_symbol=tc.get("gene_symbol", ""),
                transcript_id=tc.get("transcript_id", ""),
                consequence=",".join(tc.get("consequence_terms", [])),
                impact=tc.get("impact", ""),
                sift=tc.get("sift_prediction"),
                polyphen=tc.get("polyphen_prediction"),
            ))
    return results


# ---------------------------------------------------------------------------
# gnomAD
# ---------------------------------------------------------------------------

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
_GNOMAD_QUERY = """
query GeneConstraint($gene: String!) {
  gene(gene_symbol: $gene, reference_genome: GRCh38) {
    gnomad_constraint {
      pLI
      oe_lof_upper
    }
  }
}
"""


def _query_gnomad_constraint(gene: str) -> Optional[dict]:
    try:
        resp = requests.post(
            GNOMAD_URL,
            json={"query": _GNOMAD_QUERY, "variables": {"gene": gene}},
            timeout=15,
        )
        if resp.status_code != 200:
            return None
        data = resp.json()
        return (
            data.get("data", {})
                .get("gene", {})
                .get("gnomad_constraint")
        )
    except requests.RequestException:
        return None


def query_gene_constraints(genes: list[str]) -> list[GeneConstraint]:
    results: list[GeneConstraint] = []
    for gene in genes:
        constraint = _query_gnomad_constraint(gene)
        if constraint is None:
            # Fallback: assume unknown gene is targetable
            results.append(GeneConstraint(gene=gene, pLI=0.0, LOEUF=1.0, safe_to_target=True))
            continue

        pLI   = constraint.get("pLI")   or 0.0
        loeuf = constraint.get("oe_lof_upper") or 1.0
        # Canonical thresholds from gnomAD paper (Karczewski et al. 2020)
        safe  = pLI < 0.9 and loeuf > 0.35
        results.append(GeneConstraint(gene=gene, pLI=pLI, LOEUF=loeuf, safe_to_target=safe))
        time.sleep(0.5)   # gentle rate-limit
    return results


# ---------------------------------------------------------------------------
# ClinVar E-utilities
# ---------------------------------------------------------------------------

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_EMAIL = "knotworking-pipeline@example.com"   # NCBI requests you identify yourself


def _clinvar_significance(variant_id: str) -> Optional[str]:
    """Look up clinical significance for an rs ID or variant name in ClinVar."""
    try:
        search = requests.get(
            ESEARCH_URL,
            params={"db": "clinvar", "term": variant_id, "retmax": 1,
                    "retmode": "json", "email": NCBI_EMAIL},
            timeout=10,
        ).json()
        ids = search.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None

        summary = requests.get(
            ESUMMARY_URL,
            params={"db": "clinvar", "id": ids[0], "retmode": "json",
                    "email": NCBI_EMAIL},
            timeout=10,
        ).json()
        result = summary.get("result", {}).get(ids[0], {})
        return result.get("clinical_significance", {}).get("description")
    except requests.RequestException:
        return None


# ---------------------------------------------------------------------------
# Main Layer 3 class
# ---------------------------------------------------------------------------

class Layer3GenoPheno:
    """
    Annotates patient variants and returns dysfunctional proteins for
    downstream pathway modelling (Layer 4).
    """

    HIGH_IMPACT = {"HIGH", "MODERATE"}
    DAMAGING    = {"deleterious", "probably_damaging", "possibly_damaging"}

    def run(self, hgvs_notations: list[str], genes: list[str]) -> Layer3Result:
        result = Layer3Result()

        # 1. Variant consequence annotation via Ensembl VEP
        if hgvs_notations:
            print(f"[L3] Querying Ensembl VEP for {len(hgvs_notations)} variant(s)...")
            raw = _query_vep(hgvs_notations)
            result.consequences = parse_vep(raw)

            # Enrich with ClinVar significance
            for vc in result.consequences:
                sig = _clinvar_significance(vc.variant_id)
                if sig:
                    vc.clinical_sig = sig

        # 2. Gene constraint scores via gnomAD
        if genes:
            print(f"[L3] Querying gnomAD constraint for {genes}...")
            result.constraints = query_gene_constraints(genes)

        # 3. Derive list of likely dysfunctional proteins
        dysfunctional: set[str] = set()
        for vc in result.consequences:
            if vc.impact in self.HIGH_IMPACT:
                if vc.sift in self.DAMAGING or vc.polyphen in self.DAMAGING:
                    if vc.gene_symbol:
                        dysfunctional.add(vc.gene_symbol)
        # Also flag genes with pathogenic ClinVar entries
        for vc in result.consequences:
            if vc.clinical_sig and "pathogenic" in vc.clinical_sig.lower():
                if vc.gene_symbol:
                    dysfunctional.add(vc.gene_symbol)

        result.dysfunctional_proteins = sorted(dysfunctional)
        return result


if __name__ == "__main__":
    l3 = Layer3GenoPheno()
    out = l3.run(
        hgvs_notations=["NM_000633.3:c.118A>G"],
        genes=["BCL2", "MYC", "TP53"],
    )
    print(f"\nDysfunctional proteins: {out.dysfunctional_proteins}")
    for gc in out.constraints:
        print(f"  {gc.gene}: pLI={gc.pLI:.3f}  LOEUF={gc.LOEUF:.3f}  safe_to_target={gc.safe_to_target}")
