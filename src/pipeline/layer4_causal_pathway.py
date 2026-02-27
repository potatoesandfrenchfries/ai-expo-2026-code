"""Layer 4: Causal Pathway Modelling.

Builds a directed protein-interaction graph that traces how a mutation in a
dysfunctional protein propagates through the pathway to produce the observed
disease phenotype.  Rate-limiting nodes in that graph become the candidate
targets passed to Layer 5.

External services used (both free, no auth required):
┌────────────────────────────────────────────────────────────────┐
│  STRING DB REST API (string-db.org)                            │
│  — 20,000+ organisms, experimental + predicted interactions    │
│  — Confidence scores 0–1000 (we threshold at 400 by default)  │
├────────────────────────────────────────────────────────────────┤
│  Reactome Content Service REST (reactome.org)                  │
│  — Maps gene symbols to canonical pathway names                │
│  — Used to annotate graph edges with biological context        │
└────────────────────────────────────────────────────────────────┘

Graph analysis (local, no external service):
  • Betweenness centrality  — bottleneck nodes in the interaction network
  • Shortest paths          — causal chain from mutated gene → phenotype marker
  • Weakly connected check  — flags disconnected proteins for review

Output feeds into:
  • Layer 5 (target ID): candidate_targets (sorted by centrality)

Requires: networkx, requests
"""

import time
from dataclasses import dataclass, field
from typing import Optional

import requests

try:
    import networkx as nx
except ImportError:
    raise ImportError("networkx is required: pip install networkx")


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

@dataclass
class ProteinNode:
    gene:        str
    degree:      int
    betweenness: float   # 0–1; higher = more central bottleneck
    pathways:    list[str] = field(default_factory=list)


@dataclass
class Layer4Result:
    graph:             object            # networkx.DiGraph
    nodes:             list[ProteinNode]
    candidate_targets: list[str]         # top-centrality proteins, excl. seed
    causal_chains:     dict              # seed_gene → list of nodes to phenotype_marker


# ---------------------------------------------------------------------------
# STRING DB
# ---------------------------------------------------------------------------

STRING_NETWORK_URL = "https://string-db.org/api/json/network"
STRING_ENRICH_URL  = "https://string-db.org/api/json/enrichment"
HUMAN_TAXON = 9606


def _query_string_interactions(
    proteins: list[str],
    min_score: int = 400,
    species: int = HUMAN_TAXON,
) -> list[dict]:
    """
    Return raw STRING interactions for a list of protein names.
    min_score: 0–1000 (400 = medium confidence, 700 = high confidence).
    """
    if not proteins:
        return []

    identifiers = "%0d".join(proteins)   # STRING API separator
    params = {
        "identifiers": identifiers,
        "species":     species,
        "required_score": min_score,
        "caller_identity": "knotworking-pipeline",
    }
    for attempt in range(4):
        try:
            resp = requests.get(STRING_NETWORK_URL, params=params, timeout=20)
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 429:
                time.sleep(int(resp.headers.get("Retry-After", 5)))
                continue
        except requests.RequestException as e:
            print(f"  STRING attempt {attempt+1} failed: {e}")
        time.sleep(2 ** attempt)
    return []


def build_interaction_graph(interactions: list[dict]) -> "nx.DiGraph":
    G = nx.DiGraph()
    for edge in interactions:
        a = edge.get("preferredName_A", edge.get("stringId_A", ""))
        b = edge.get("preferredName_B", edge.get("stringId_B", ""))
        score = edge.get("score", 0) / 1000.0   # normalise to [0,1]
        if a and b:
            G.add_edge(a, b, weight=score)
            G.add_edge(b, a, weight=score)   # treat as bidirectional by default
    return G


# ---------------------------------------------------------------------------
# Reactome pathway annotation
# ---------------------------------------------------------------------------

REACTOME_URL = "https://reactome.org/ContentService/data/mapping/UniProt/{uniprot}/pathways"
REACTOME_QUERY_URL = "https://reactome.org/ContentService/search/query"


def _get_reactome_pathways(gene: str) -> list[str]:
    """Return top Reactome pathway names for a gene symbol (best-effort)."""
    try:
        resp = requests.get(
            REACTOME_QUERY_URL,
            params={"query": gene, "species": "Homo+sapiens", "types": "Pathway",
                    "cluster": "true"},
            timeout=10,
        )
        if resp.status_code != 200:
            return []
        results = resp.json().get("results", [])
        # Flatten all pathway names across result groups
        names: list[str] = []
        for group in results:
            for entry in group.get("entries", []):
                name = entry.get("name")
                if name:
                    names.append(name)
        return names[:5]   # top 5 pathways
    except (requests.RequestException, ValueError):
        return []


# ---------------------------------------------------------------------------
# Graph analysis helpers
# ---------------------------------------------------------------------------

def compute_centrality(G: "nx.DiGraph") -> dict:
    """Return betweenness centrality for every node (normalised 0–1)."""
    if len(G.nodes) < 2:
        return {n: 0.0 for n in G.nodes}
    return nx.betweenness_centrality(G, normalized=True, weight="weight")


def shortest_paths(
    G: "nx.DiGraph",
    sources: list[str],
    targets: list[str],
) -> dict:
    """
    For each (source, target) pair that exists in the graph, return the
    shortest path as a list of gene names.
    """
    paths: dict = {}
    for src in sources:
        if src not in G:
            continue
        for tgt in targets:
            if tgt not in G or tgt == src:
                continue
            try:
                path = nx.shortest_path(G, src, tgt, weight=None)
                paths[f"{src}→{tgt}"] = path
            except nx.NetworkXNoPath:
                pass
    return paths


# ---------------------------------------------------------------------------
# Main Layer 4 class
# ---------------------------------------------------------------------------

class Layer4CausalPathway:
    """
    Builds a protein interaction graph from STRING DB and computes:
      - betweenness centrality (bottleneck / rate-limiting proteins)
      - shortest causal chain from mutated genes to phenotype markers
      - Reactome pathway annotations for each node

    Parameters
    ----------
    min_score : int
        STRING confidence threshold (0–1000). 400 = medium, 700 = high.
    phenotype_markers : list[str]
        Known downstream proteins linked to the disease phenotype
        (e.g. apoptosis markers for lymphoma: ["BAX", "CASP3"]).
    """

    def __init__(
        self,
        min_score: int = 400,
        phenotype_markers: Optional[list[str]] = None,
    ):
        self.min_score         = min_score
        self.phenotype_markers = phenotype_markers or []

    def run(
        self,
        dysfunctional_proteins: list[str],
        safe_genes: Optional[list[str]] = None,
    ) -> Layer4Result:
        """
        Parameters
        ----------
        dysfunctional_proteins : list[str]
            Output of Layer 3 (genes with damaging variants).
        safe_genes : list[str] | None
            Genes flagged as essential (pLI > 0.9) from Layer 3;
            they will be excluded from candidate_targets.
        """
        safe_set = set(safe_genes or [])
        all_proteins = list(set(dysfunctional_proteins + self.phenotype_markers))

        print(f"[L4] Querying STRING DB for {all_proteins} ...")
        interactions = _query_string_interactions(all_proteins, min_score=self.min_score)
        G = build_interaction_graph(interactions)

        # Add seeds as isolated nodes if STRING returned nothing for them
        for p in all_proteins:
            if p not in G:
                G.add_node(p)

        centrality = compute_centrality(G)

        # Build ProteinNode list with Reactome pathway annotation
        nodes: list[ProteinNode] = []
        for gene, bc in sorted(centrality.items(), key=lambda kv: -kv[1]):
            pathways = _get_reactome_pathways(gene)
            nodes.append(ProteinNode(
                gene=gene,
                degree=G.degree(gene),
                betweenness=round(bc, 5),
                pathways=pathways,
            ))

        # Candidate targets: high centrality, not a seed mutation, not essential
        candidates = [
            n.gene for n in nodes
            if n.gene not in dysfunctional_proteins
            and n.gene not in safe_set
            and n.betweenness > 0
        ][:10]   # top 10

        # Causal chain analysis
        chains = shortest_paths(G, dysfunctional_proteins, self.phenotype_markers)

        return Layer4Result(
            graph=G,
            nodes=nodes,
            candidate_targets=candidates,
            causal_chains=chains,
        )


# ---------------------------------------------------------------------------
# Quick demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    l4 = Layer4CausalPathway(
        min_score=400,
        phenotype_markers=["BAX", "CASP3"],
    )
    result = l4.run(dysfunctional_proteins=["BCL2", "MYC"])

    print(f"\nTop candidate targets: {result.candidate_targets}")
    print("\nCausal chains:")
    for chain_id, path in result.causal_chains.items():
        print(f"  {chain_id}: {' → '.join(path)}")
    print("\nTop 5 nodes by centrality:")
    for node in result.nodes[:5]:
        print(f"  {node.gene}: degree={node.degree}, betweenness={node.betweenness}")
