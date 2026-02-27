"""GitHub Rocq/Coq corpus scraper.

Searches GitHub for .v files containing chemistry and drug-discovery proof
patterns and downloads them to data/github_rocq_corpus/.  The scraped files
are then used by training_data_generator.py to expand the few-shot proof
corpus beyond the 13 hand-written Chemistry library modules.

Authentication:
    Set GITHUB_TOKEN env var for 30 req/min (vs 10 req/min unauthenticated).

Usage:
    python -m src.data.github_scraper [--max-files 200] [--out-dir data/github_rocq_corpus]
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

# Rocq/Coq proof patterns relevant to chemistry + drug-likeness
SEARCH_QUERIES = [
    "language:coq molecule drug extension:v",
    "language:coq Lipinski pharmacology extension:v",
    "language:coq atom bond valence extension:v",
    "language:coq chemistry proof extension:v",
    "language:coq SMILES molecular extension:v",
]

# Must contain at least one of these to be considered chemistry-relevant
RELEVANCE_PATTERNS = [
    r"\bMolecule\b", r"\bAtom\b", r"\bBond\b",
    r"\bLipinski\b", r"\bvalence\b", r"\bdrug\b",
    r"\bpharma\b", r"\bSMILES\b", r"\bchemistry\b",
]

_RELEVANCE_RE = re.compile("|".join(RELEVANCE_PATTERNS), re.IGNORECASE)


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


def download_file(item: dict, out_dir: Path) -> Path | None:
    """Download a single .v file and return its local path, or None on failure."""
    repo     = item["repository"]["full_name"]
    path     = item["path"]
    sha      = item["sha"]
    ref      = item.get("ref", "HEAD")

    # Sanitise filename: replace directory separators with __
    safe_name = f"{repo.replace('/', '__')}__{path.replace('/', '__')}"
    out_path  = out_dir / safe_name

    if out_path.exists():
        return out_path  # already cached

    raw_url = f"{GITHUB_RAW_URL}/{repo}/{ref}/{path}"
    for attempt in range(4):
        resp = requests.get(raw_url, headers=_headers(), timeout=15)
        if resp.status_code == 200:
            break
        _rate_limit_wait(resp, attempt)
    else:
        return None

    content = resp.text
    if not _RELEVANCE_RE.search(content):
        return None  # not chemistry-relevant, skip

    out_path.write_text(content, encoding="utf-8")
    return out_path


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


def extract_proof_patterns(v_file: Path) -> list[dict]:
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
        patterns.append({"type": "lemma", "source": str(v_file.name), "content": chunk})

    for m in _DEFINITION_RE.finditer(text):
        chunk = m.group(1).strip()
        if len(chunk) < 20 or len(chunk) > 800:
            continue
        patterns.append({"type": "definition", "source": str(v_file.name), "content": chunk})

    return patterns


# ---------------------------------------------------------------------------
# Main entrypoint
# ---------------------------------------------------------------------------

def scrape(max_files: int = 150, out_dir: str = "data/github_rocq_corpus") -> None:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    patterns_path = Path("data") / "github_proof_patterns.jsonl"

    seen_shas: set[str] = set()
    downloaded: list[Path] = []
    all_patterns: list[dict] = []

    files_per_query = max(1, max_files // len(SEARCH_QUERIES))

    for query in SEARCH_QUERIES:
        print(f"\nSearching: '{query}'")
        items = search_rocq_files(query, max_results=files_per_query)
        print(f"  Found {len(items)} candidate files.")

        for item in items:
            sha = item["sha"]
            if sha in seen_shas:
                continue
            seen_shas.add(sha)

            local = download_file(item, out_path)
            if local is None:
                continue

            downloaded.append(local)
            patterns = extract_proof_patterns(local)
            all_patterns.extend(patterns)
            print(f"  + {local.name}  ({len(patterns)} patterns)")

            time.sleep(0.5)  # gentle secondary rate-limit buffer

    # Write extracted patterns to JSONL
    with open(patterns_path, "w") as f:
        for p in all_patterns:
            f.write(json.dumps(p) + "\n")

    print(f"\nDone. {len(downloaded)} files, {len(all_patterns)} proof patterns.")
    print(f"Files    → {out_path}/")
    print(f"Patterns → {patterns_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scrape Rocq chemistry proofs from GitHub.")
    parser.add_argument("--max-files", type=int, default=150)
    parser.add_argument("--out-dir",   type=str,  default="data/github_rocq_corpus")
    args = parser.parse_args()
    scrape(max_files=args.max_files, out_dir=args.out_dir)
