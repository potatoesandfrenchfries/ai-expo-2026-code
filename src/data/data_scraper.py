"""GitHub Rocq/Coq corpus scraper.

Searches GitHub for .v files across three domains and downloads them to
data/github_rocq_corpus/.  The scraped files are used by
training_data_generator.py to expand the few-shot proof corpus.

Domain strategy
---------------
Very few public repos formalise chemistry directly in Rocq/Coq.  However,
our Chemistry library proofs are entirely built on:
  • Real-number arithmetic  (lra, field_simplify, Reals library)
  • Linear algebra / vectors  (used in Geometry.v)
  • Graph theory  (used in Molecules.v, MathProperties.v)
  • Inequalities and bounds  (used throughout DrugLikeness.v, Valency.v)
  • Probability / measure  (Bayesian uncertainty layer)
  • Physics thermodynamics  (energy, entropy — relevant to binding affinity)

All of these have substantial Rocq corpora on GitHub.  Scraping them gives
the proof LLM a much richer diet of tactic patterns (lra, lia, ring, auto,
omega, field, norm_num) and proof structures to learn from.

Authentication:
    Set GITHUB_TOKEN env var for 30 req/min (vs 10 req/min unauthenticated).

Usage:
    python -m src.data.github_scraper [--max-files 300] [--out-dir data/github_rocq_corpus]
"""

import argparse
import os
import time
import json
import re
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

GITHUB_SEARCH_URL = "https://api.github.com/search/code"
GITHUB_RAW_URL    = "https://raw.githubusercontent.com"
GITHUB_CONTENTS_URL = "https://api.github.com/repos/{repo}/git/trees/{sha}?recursive=1"

# ── Search queries grouped by domain ────────────────────────────────────────
# Each entry: (query_string, domain_tag)
SEARCH_QUERIES: list[tuple[str, str]] = [
    # Chemistry / drug-discovery (direct)
    ("language:coq molecule drug extension:v",               "chemistry"),
    ("language:coq Lipinski pharmacology extension:v",       "chemistry"),
    ("language:coq atom bond valence extension:v",           "chemistry"),
    ("language:coq SMILES molecular extension:v",            "chemistry"),
    # Real analysis — foundation of our numeric bounds proofs
    ("language:coq real analysis lra extension:v",           "math_analysis"),
    ("language:coq Reals field_simplify extension:v",        "math_analysis"),
    ("language:coq norm_num inequality extension:v",         "math_analysis"),
    ("language:coq integral differential extension:v",       "math_analysis"),
    # Linear algebra / vectors — Geometry.v, rotation matrices
    ("language:coq vector matrix linear_map extension:v",    "math_linalg"),
    ("language:coq inner_product norm metric extension:v",   "math_linalg"),
    # Graph theory — Molecules.v connectivity, BFS
    ("language:coq graph path connected extension:v",        "math_graph"),
    ("language:coq cycle subgraph isomorphism extension:v",  "math_graph"),
    # Probability / measure — Bayesian uncertainty
    ("language:coq probability measure sigma extension:v",   "math_prob"),
    ("language:coq distribution expectation variance extension:v", "math_prob"),
    # Physics / thermodynamics — energy, entropy relevant to binding
    ("language:coq energy entropy thermodynamics extension:v", "physics"),
    ("language:coq quantum Hamiltonian extension:v",         "physics"),
    ("language:coq oscillator potential force extension:v",  "physics"),
]

# ── Known high-value repos to scrape directly (bypass search limits) ────────
# These are the canonical Coq/Rocq maths libraries whose tactic patterns
# appear throughout our Chemistry proofs.
DIRECT_REPOS: list[tuple[str, str, str]] = [
    # (owner/repo, branch, domain_tag)
    ("math-comp/math-comp",          "master",  "math_analysis"),   # algebra, real analysis
    ("math-comp/analysis",           "master",  "math_analysis"),   # Coquelicot-style analysis
    ("coq-community/graph-theory",   "master",  "math_graph"),      # graph theory
    ("coq-community/topology",       "master",  "math_analysis"),   # metric spaces
    ("coq-community/coq-ext-lib",    "master",  "math_analysis"),   # utilities
]

# ── Relevance filter ─────────────────────────────────────────────────────────
# A file must match at least one pattern from its domain group to be kept.
_DOMAIN_PATTERNS: dict[str, list[str]] = {
    "chemistry": [
        r"\bMolecule\b", r"\bAtom\b", r"\bBond\b", r"\bLipinski\b",
        r"\bvalence\b", r"\bdrug\b", r"\bpharma\b", r"\bSMILES\b",
    ],
    "math_analysis": [
        r"\bRbar\b", r"\bRlist\b", r"\blra\b", r"\bfield_simplify\b",
        r"\bnorm_num\b", r"\bIntegral\b", r"\bDerivative\b",
        r"\bContinuous\b", r"\bConverge\b", r"\bMetric\b",
        r"\bR\b.*\bProof\b",   # proofs about real numbers
    ],
    "math_linalg": [
        r"\bMatrix\b", r"\bVector\b", r"\binner_product\b",
        r"\blinear_map\b", r"\beigenvalue\b", r"\borthogonal\b",
        r"\bnorm\b", r"\bdot_product\b",
    ],
    "math_graph": [
        r"\bGraph\b", r"\bPath\b", r"\bConnected\b", r"\bCycle\b",
        r"\bBFS\b", r"\bDFS\b", r"\bSubgraph\b", r"\bIsomorphism\b",
        r"\bDegree\b.*\bvertex\b",
    ],
    "math_prob": [
        r"\bProbability\b", r"\bMeasure\b", r"\bSigma\b",
        r"\bDistribution\b", r"\bExpectation\b", r"\bVariance\b",
        r"\bBayes\b", r"\bConditional\b",
    ],
    "physics": [
        r"\bEnergy\b", r"\bEntropy\b", r"\bHamiltonian\b",
        r"\bthermodynamic\b", r"\bQuantum\b", r"\boscillator\b",
        r"\bpotential\b", r"\bforce\b", r"\btemperature\b",
    ],
}

# Fallback: any file with a Proof block is weakly relevant
_FALLBACK_RE = re.compile(r"\bProof\b.*\bQed\b", re.DOTALL)

def _is_relevant(content: str, domain: str) -> bool:
    patterns = _DOMAIN_PATTERNS.get(domain, [])
    combined = re.compile("|".join(patterns), re.IGNORECASE) if patterns else None
    if combined and combined.search(content):
        return True
    return bool(_FALLBACK_RE.search(content))


# ---------------------------------------------------------------------------
# API helpers
# ---------------------------------------------------------------------------

def _headers() -> dict:
    token = os.environ.get("GITHUB_TOKEN")
    h = {"Accept": "application/vnd.github+json"}
    if token:
        h["Authorization"] = f"Bearer {token}"
    return h


def _rate_limit_wait(response: requests.Response, attempt: int) -> None:
    """Respect GitHub rate limits and Retry-After headers."""
    if response.status_code == 403:
        reset = int(response.headers.get("X-RateLimit-Reset", time.time() + 60))
        wait = max(reset - int(time.time()), 1)
        print(f"  Rate limited. Waiting {wait}s ...")
        time.sleep(wait)
    else:
        # Exponential backoff for other transient errors
        time.sleep(min(2 ** attempt, 30))


def search_rocq_files(query: str, max_results: int = 100) -> list[dict]:
    """Return a list of GitHub file metadata dicts for a given search query."""
    items: list[dict] = []
    page = 1
    per_page = 30  # GitHub caps code search at 30/page

    while len(items) < max_results:
        params = {"q": query, "per_page": per_page, "page": page}
        for attempt in range(4):
            resp = requests.get(GITHUB_SEARCH_URL, headers=_headers(), params=params, timeout=15)
            if resp.status_code == 200:
                break
            print(f"  Search attempt {attempt+1} failed ({resp.status_code}): {resp.text[:80]}")
            _rate_limit_wait(resp, attempt)
        else:
            print(f"  Giving up on query page {page}.")
            break

        data = resp.json()
        batch = data.get("items", [])
        if not batch:
            break

        items.extend(batch)
        page += 1

        # Respect secondary rate limits between pages
        time.sleep(1.0)

        if len(items) >= data.get("total_count", 0):
            break

    return items[:max_results]


def download_file(item: dict, out_dir: Path, domain: str = "chemistry") -> tuple:
    """
    Download a single .v file, apply domain relevance filter.
    Returns (local_path, domain) or (None, domain) on failure/irrelevance.
    """
    repo  = item["repository"]["full_name"]
    path  = item["path"]
    ref   = item.get("ref", "HEAD")

    safe_name = f"{repo.replace('/', '__')}__{path.replace('/', '__')}"
    out_path  = out_dir / safe_name

    if out_path.exists():
        return out_path, domain

    raw_url = f"{GITHUB_RAW_URL}/{repo}/{ref}/{path}"
    for attempt in range(4):
        resp = requests.get(raw_url, headers=_headers(), timeout=15)
        if resp.status_code == 200:
            break
        _rate_limit_wait(resp, attempt)
    else:
        return None, domain

    content = resp.text
    if not _is_relevant(content, domain):
        return None, domain

    out_path.write_text(content, encoding="utf-8")
    return out_path, domain


def scrape_repo_directly(
    repo: str, branch: str, domain: str, out_dir: Path, max_files: int = 30
) -> list[tuple]:
    """
    Walk a known high-value repo's file tree and download all .v files
    up to max_files.  Bypasses search API so it doesn't count against
    the 10/30 req/min code-search quota.
    """
    # Get the tree SHA for the branch
    branch_url = f"https://api.github.com/repos/{repo}/branches/{branch}"
    resp = requests.get(branch_url, headers=_headers(), timeout=15)
    if resp.status_code != 200:
        print(f"  Could not fetch branch info for {repo}: {resp.status_code}")
        return []

    tree_sha = resp.json().get("commit", {}).get("commit", {}).get("tree", {}).get("sha")
    if not tree_sha:
        return []

    # Recursive tree listing
    tree_url = f"https://api.github.com/repos/{repo}/git/trees/{tree_sha}?recursive=1"
    resp = requests.get(tree_url, headers=_headers(), timeout=20)
    if resp.status_code != 200:
        return []

    v_files = [
        item for item in resp.json().get("tree", [])
        if item.get("path", "").endswith(".v") and item.get("type") == "blob"
    ][:max_files]

    results = []
    for tree_item in v_files:
        fake_item = {
            "repository": {"full_name": repo},
            "path": tree_item["path"],
            "sha":  tree_item["sha"],
            "ref":  branch,
        }
        local, dom = download_file(fake_item, out_dir, domain)
        if local:
            results.append((local, dom))
        time.sleep(0.3)

    return results


# ---------------------------------------------------------------------------
# Proof pattern extractor
# ---------------------------------------------------------------------------

_LEMMA_RE = re.compile(
    r"((?:Lemma|Theorem|Proposition|Corollary)\s+\w+.*?(?:Qed|Defined|Admitted)\.)",
    re.DOTALL,
)
_DEFINITION_RE = re.compile(
    r"(Definition\s+\w+.*?\.)",
    re.DOTALL,
)


def extract_proof_patterns(v_file: Path, domain: str = "unknown") -> list[dict]:
    """
    Extract (statement, proof) pairs from a .v file.
    Each pair becomes a candidate training example for the proof LLM.
    """
    text = v_file.read_text(encoding="utf-8", errors="replace")
    patterns: list[dict] = []

    for m in _LEMMA_RE.finditer(text):
        chunk = m.group(1).strip()
        if len(chunk) < 30 or len(chunk) > 2000:
            continue
        patterns.append({"type": "lemma", "domain": domain, "source": str(v_file.name), "content": chunk})

    for m in _DEFINITION_RE.finditer(text):
        chunk = m.group(1).strip()
        if len(chunk) < 20 or len(chunk) > 800:
            continue
        patterns.append({"type": "definition", "domain": domain, "source": str(v_file.name), "content": chunk})

    return patterns


# ---------------------------------------------------------------------------
# Main entrypoint
# ---------------------------------------------------------------------------

def scrape(max_files: int = 300, out_dir: str = "data/github_rocq_corpus") -> None:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    patterns_path = Path("data") / "github_proof_patterns.jsonl"

    seen_shas: set[str] = set()
    downloaded: list[tuple] = []   # (path, domain)
    all_patterns: list[dict] = []

    # ── Phase 1: Direct scrape of known high-value repos ────────────────────
    print("=== Phase 1: Direct repo scraping ===")
    files_per_repo = max(10, max_files // (len(DIRECT_REPOS) * 4))
    for repo, branch, domain in DIRECT_REPOS:
        print(f"\nScraping {repo} ({domain}) ...")
        results = scrape_repo_directly(repo, branch, domain, out_path, max_files=files_per_repo)
        for local, dom in results:
            if local not in [d[0] for d in downloaded]:
                downloaded.append((local, dom))
                patterns = extract_proof_patterns(local, domain=dom)
                all_patterns.extend(patterns)
                print(f"  + {local.name}  ({len(patterns)} patterns)")
        time.sleep(1.0)

    # ── Phase 2: Search API across all domain queries ────────────────────────
    print("\n=== Phase 2: GitHub code search ===")
    files_per_query = max(1, (max_files - len(downloaded)) // len(SEARCH_QUERIES))

    for query, domain in SEARCH_QUERIES:
        print(f"\nSearching [{domain}]: '{query}'")
        items = search_rocq_files(query, max_results=files_per_query)
        print(f"  Found {len(items)} candidate files.")

        for item in items:
            sha = item["sha"]
            if sha in seen_shas:
                continue
            seen_shas.add(sha)

            local, dom = download_file(item, out_path, domain)
            if local is None:
                continue

            downloaded.append((local, dom))
            patterns = extract_proof_patterns(local, domain=dom)
            all_patterns.extend(patterns)
            print(f"  + {local.name}  ({len(patterns)} patterns)")

            time.sleep(0.5)

    # ── Write corpus ─────────────────────────────────────────────────────────
    with open(patterns_path, "w") as f:
        for p in all_patterns:
            f.write(json.dumps(p) + "\n")

    # Domain breakdown summary
    domain_counts: dict[str, int] = {}
    for p in all_patterns:
        d = p.get("domain", "unknown")
        domain_counts[d] = domain_counts.get(d, 0) + 1

    print(f"\nDone. {len(downloaded)} files, {len(all_patterns)} proof patterns.")
    print("Breakdown by domain:")
    for dom, cnt in sorted(domain_counts.items(), key=lambda x: -x[1]):
        print(f"  {dom:<20} {cnt:>5} patterns")
    print(f"\nFiles    → {out_path}/")
    print(f"Patterns → {patterns_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Scrape Rocq proofs from GitHub: chemistry, maths, and physics."
    )
    parser.add_argument("--max-files", type=int, default=300)
    parser.add_argument("--out-dir",   type=str,  default="data/github_rocq_corpus")
    args = parser.parse_args()
    scrape(max_files=args.max_files, out_dir=args.out_dir)