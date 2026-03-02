"""Rocq subprocess bridge — structured per-constraint proof checking.

Provides check_molecule(smiles) → dict, a pure-function interface that:
  1. Checks PAINS / SA score quickly via RDKit (no subprocess).
  2. Generates a Rocq proof file encoding all numeric constraints
     (Lipinski RO5 + Veber rules) as individual named theorems.
  3. Runs coqc (or coqtop -batch) to formally compile the proof.
  4. Parses the subprocess output to identify the specific failing
     lemma name, its hypothesis string, and the violation magnitude.
  5. Returns a structured result dict; caches results by InChIKey.

The bridge treats the Rocq proof checker as a black-box oracle — no
gradients flow through it.  It is called from constraint_vector.py and
rl_reward.py.

Reference: Lipinski et al. 2001; Baell & Holloway 2010 (J Med Chem 53:2719).
"""

import os
import re
import subprocess
import tempfile
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# Optional SA score (rdkit.Contrib may not be installed in all environments)
try:
    from rdkit.Contrib.SA_Score import sascorer as _sascorer
    _HAS_SA_SCORE = True
except Exception:
    _HAS_SA_SCORE = False

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)
_ROCQ_LIB = os.path.join(_REPO_ROOT, "src", "rocq")

# ---------------------------------------------------------------------------
# In-process cache: InChIKey → result dict
# ---------------------------------------------------------------------------

_proof_cache: dict = {}


def check_molecule(smiles: str, timeout: float = 10.0) -> dict:
    """Check a SMILES string against all drug-likeness constraints.

    Returns:
        {
            "passed": bool,
            "failing_lemma": str | None,       # e.g. "lipinski_logp_constraint"
            "failing_hypothesis": str | None,  # e.g. "logp_constraint : 6.30 > 5"
            "constraint_type": str | None,     # one of: mw logp hbd hba
                                               #         pains sa_score toxicity
                                               #         (rot psa also possible)
            "violation_magnitude": float | None,  # how far outside the bound
        }

    The call is idempotent and cached by InChIKey so redundant Rocq
    subprocess invocations are avoided.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "passed": False,
            "failing_lemma": "invalid_smiles",
            "failing_hypothesis": f"SMILES could not be parsed: {smiles!r}",
            "constraint_type": None,
            "violation_magnitude": None,
        }

    # --- cache lookup ---
    inchi_key = Chem.MolToInchiKey(mol) or ""
    if inchi_key and inchi_key in _proof_cache:
        return _proof_cache[inchi_key]

    # --- fast RDKit filters (no subprocess) ---
    pains_result = _check_pains(mol)
    if pains_result is not None:
        if inchi_key:
            _proof_cache[inchi_key] = pains_result
        return pains_result

    sa_result = _check_sa_score(mol)
    if sa_result is not None:
        if inchi_key:
            _proof_cache[inchi_key] = sa_result
        return sa_result

    # --- formal Rocq verification of numeric constraints ---
    desc = _compute_descriptors(mol)
    result = _run_rocq_check(desc, timeout)

    if inchi_key:
        _proof_cache[inchi_key] = result
    return result


def clear_cache() -> None:
    """Discard all cached proof results (e.g. between training epochs)."""
    _proof_cache.clear()


def cache_size() -> int:
    """Return the number of cached InChIKey entries."""
    return len(_proof_cache)


# ---------------------------------------------------------------------------
# RDKit fast checks
# ---------------------------------------------------------------------------

def _compute_descriptors(mol) -> dict:
    logp = Descriptors.MolLogP(mol)
    return {
        "mw":            round(Descriptors.ExactMolWt(mol), 4),
        "logp":          round(logp, 4),
        "hbd":           rdMolDescriptors.CalcNumHBD(mol),
        "hba":           rdMolDescriptors.CalcNumHBA(mol),
        "rot":           rdMolDescriptors.CalcNumRotatableBonds(mol),
        "psa":           round(Descriptors.TPSA(mol), 4),
        "molar_refract": round(Descriptors.MolMR(mol), 4),
        "logS":          round(-logp * 0.54 - 0.54, 4),
        "rings":         rdMolDescriptors.CalcNumRings(mol),
        "arom_rings":    rdMolDescriptors.CalcNumAromaticRings(mol),
        "chiral":        len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
        "heavy_atoms":   mol.GetNumHeavyAtoms(),
    }


def _check_pains(mol) -> Optional[dict]:
    """Return a failure dict if a PAINS or Brenk alert matches, else None."""
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog(params)
    entry = catalog.GetFirstMatch(mol)
    if entry is None:
        return None
    alert_name = entry.GetDescription()
    safe_name = re.sub(r"[^a-z0-9]+", "_", alert_name.lower()).strip("_")
    return {
        "passed": False,
        "failing_lemma": f"pains_{safe_name}_filter",
        "failing_hypothesis": f"passes_pains : alert={alert_name!r}",
        "constraint_type": "pains",
        "violation_magnitude": 1.0,
    }


def _check_sa_score(mol) -> Optional[dict]:
    """Return a failure dict if SA score > 6.0 (hard to synthesise), else None."""
    if not _HAS_SA_SCORE:
        return None
    try:
        sa = _sascorer.calculateScore(mol)
        if sa > 6.0:
            return {
                "passed": False,
                "failing_lemma": "sa_score_constraint",
                "failing_hypothesis": f"sa_score_under_6 : {sa:.3f} > 6.0",
                "constraint_type": "sa_score",
                "violation_magnitude": round(sa - 6.0, 4),
            }
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Rocq proof generation
# ---------------------------------------------------------------------------

# Ordered list of numeric constraints to verify formally.
# Each entry: (constraint_type, theorem_name, Rocq proposition, tactic)
_NUMERIC_CONSTRAINTS = [
    ("mw",   "lipinski_mw_constraint",
     "mol_desc.(md_mol_weight) <= 500",       "lra"),
    ("logp", "lipinski_logp_constraint",
     "mol_desc.(md_logP) <= 5",               "lra"),
    ("hbd",  "lipinski_hbd_constraint",
     "(mol_desc.(md_hbd) <= 5)%nat",          "simpl. lia"),
    ("hba",  "lipinski_hba_constraint",
     "(mol_desc.(md_hba) <= 10)%nat",         "simpl. lia"),
    ("rot",  "veber_rot_bonds_constraint",
     "(mol_desc.(md_rot_bonds) <= 10)%nat",   "simpl. lia"),
    ("psa",  "veber_psa_constraint",
     "mol_desc.(md_psa) <= 140",              "lra"),
]

_HEADER = """\
From Coq Require Import List.
Import ListNotations.
Require Import Coq.Reals.Reals.
Open Scope R_scope.
Require Import Coq.ZArith.ZArith.
Require Import Chemistry.DrugLikeness.
Require Import Coq.micromega.Lra.
Require Import Coq.micromega.Lia.
"""


def _rocq_real(v: float) -> str:
    """Format a Python float for Rocq: negative values need parentheses."""
    s = f"{v}"
    return f"({s})" if v < 0 else s


def _generate_proof(desc: dict) -> tuple[str, dict]:
    """Build the Rocq source and a mapping: tactic_line_number → (ctype, tname).

    Returns (proof_text, line_map) where line_map keys are 1-indexed line
    numbers of the tactic call inside each theorem's Proof block.  coqc
    reports the location of the failing tactic in its error message.
    """
    lines = _HEADER.splitlines()
    lines.append("")

    # mkMolDesc field order matches DrugLikeness.v record definition:
    #   mw logP hbd hba rot_bonds psa molar_refract logS rings arom_rings chiral atom_count
    moldef_val = (
        f"  mkMolDesc"
        f" {_rocq_real(desc['mw'])} {_rocq_real(desc['logp'])}"
        f" {desc['hbd']} {desc['hba']}"
        f" {desc['rot']} {_rocq_real(desc['psa'])}"
        f" {_rocq_real(desc['molar_refract'])} {_rocq_real(desc['logS'])}"
        f" {desc['rings']} {desc['arom_rings']}"
        f" {desc['chiral']} {desc['heavy_atoms']}."
    )
    lines.append("Definition mol_desc : MolDescriptors :=")
    lines.append(moldef_val)
    lines.append("")

    line_map: dict = {}  # tactic_line → (ctype, tname)

    for ctype, tname, prop, tactic in _NUMERIC_CONSTRAINTS:
        lines.append(f"(* Constraint: {ctype} *)")
        lines.append(f"Theorem {tname} : {prop}.")
        lines.append("Proof.")
        # Record 1-indexed tactic line (len(lines) counts 0-indexed, +1 for 1-indexing)
        tactic_lineno = len(lines) + 1
        lines.append(f"  unfold mol_desc. simpl. {tactic}.")
        line_map[tactic_lineno] = (ctype, tname)
        # The Qed line is at tactic_lineno + 1; record it too for robustness
        line_map[tactic_lineno + 1] = (ctype, tname)
        lines.append("Qed.")
        lines.append("")

    return "\n".join(lines) + "\n", line_map


# ---------------------------------------------------------------------------
# Rocq subprocess invocation
# ---------------------------------------------------------------------------

def _run_rocq_check(desc: dict, timeout: float) -> dict:
    """Compile the constraint proof, parse failures, and return a result dict."""
    proof_text, line_map = _generate_proof(desc)

    # Write the temp file into _ROCQ_LIB so `-R _ROCQ_LIB Chemistry` resolves
    # imports correctly (same pattern as rocq_verifier.py).
    os.makedirs(_ROCQ_LIB, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(suffix=".v", dir=_ROCQ_LIB)
    try:
        with os.fdopen(fd, "w") as f:
            f.write(proof_text)

        result = _invoke_rocq(tmp_path, timeout)

        if result is None:
            # Rocq not installed; fall back to pure Python violation check
            return _python_fallback(desc)

        if result.returncode == 0:
            return _passed()

        combined = (result.stderr or "") + (result.stdout or "")
        ctype, tname, hyp = _parse_failing(combined, line_map, desc)
        magnitude = _violation_magnitude(ctype, desc) if ctype else None
        return {
            "passed": False,
            "failing_lemma": tname,
            "failing_hypothesis": hyp,
            "constraint_type": ctype,
            "violation_magnitude": magnitude,
        }

    except subprocess.TimeoutExpired:
        return {
            "passed": False,
            "failing_lemma": "rocq_timeout",
            "failing_hypothesis": f"proof timed out after {timeout:.0f}s",
            "constraint_type": None,
            "violation_magnitude": None,
        }
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def _invoke_rocq(proof_path: str, timeout: float):
    """Try coqtop -batch, then coqc. Returns CompletedProcess or None if unavailable."""
    base_flags = ["-R", _ROCQ_LIB, "Chemistry"]

    for cmd in [
        ["coqtop", "-batch"] + base_flags + ["-l", proof_path],
        ["coqc"] + base_flags + [proof_path],
        ["rocq", "compile"] + base_flags + [proof_path],
    ]:
        try:
            return subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except FileNotFoundError:
            continue
        except subprocess.TimeoutExpired:
            raise

    return None  # no Rocq executable found


# ---------------------------------------------------------------------------
# Failure parsing
# ---------------------------------------------------------------------------

def _parse_failing(
    output: str,
    line_map: dict,
    desc: dict,
) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """Identify which constraint failed from coqc / coqtop output.

    Strategy 1 — scan for theorem name strings in the error output.
    Strategy 2 — parse the line number from the error and look it up.
    Strategy 3 — infer from descriptor values directly (always succeeds).
    """
    # Strategy 1: theorem name in output
    for ctype, tname, _, _ in _NUMERIC_CONSTRAINTS:
        if tname in output:
            hyp = _build_hypothesis(ctype, desc)
            return ctype, tname, hyp

    # Strategy 2: line number mapping
    # coqc format: "file.v:LINE:COL-COL:" or "line LINE"
    for pat in [r":(\d+):\d+[-:]", r":(\d+):", r"line (\d+)"]:
        m = re.search(pat, output)
        if m:
            error_line = int(m.group(1))
            # Walk line_map in reverse to find the nearest enclosing theorem
            for lineno in sorted(line_map.keys(), reverse=True):
                if lineno <= error_line:
                    ctype, tname = line_map[lineno]
                    hyp = _build_hypothesis(ctype, desc)
                    return ctype, tname, hyp

    # Strategy 3: pure Python fallback
    violations = _python_violations(desc)
    if violations:
        ctype, tname = violations[0]
        hyp = _build_hypothesis(ctype, desc)
        return ctype, tname, hyp

    return None, None, None


def _python_violations(desc: dict) -> list:
    """Return (ctype, tname) pairs for all violated numeric constraints."""
    checks = [
        ("mw",   "lipinski_mw_constraint",    desc["mw"]   > 500),
        ("logp", "lipinski_logp_constraint",   desc["logp"] > 5),
        ("hbd",  "lipinski_hbd_constraint",    desc["hbd"]  > 5),
        ("hba",  "lipinski_hba_constraint",    desc["hba"]  > 10),
        ("rot",  "veber_rot_bonds_constraint", desc["rot"]  > 10),
        ("psa",  "veber_psa_constraint",       desc["psa"]  > 140),
    ]
    return [(c, t) for c, t, violated in checks if violated]


def _python_fallback(desc: dict) -> dict:
    """Return a result dict derived from pure Python checks (no Rocq)."""
    violations = _python_violations(desc)
    if not violations:
        return _passed()
    ctype, tname = violations[0]
    return {
        "passed": False,
        "failing_lemma": tname,
        "failing_hypothesis": _build_hypothesis(ctype, desc),
        "constraint_type": ctype,
        "violation_magnitude": _violation_magnitude(ctype, desc),
    }


def _build_hypothesis(ctype: str, desc: dict) -> str:
    bounds = {"mw": 500, "logp": 5, "hbd": 5, "hba": 10, "rot": 10, "psa": 140}
    vals = {
        "mw": desc["mw"], "logp": desc["logp"], "hbd": desc["hbd"],
        "hba": desc["hba"], "rot": desc["rot"], "psa": desc["psa"],
    }
    if ctype in bounds:
        return f"{ctype}_constraint : {vals[ctype]} > {bounds[ctype]}"
    return f"{ctype}_constraint"


def _violation_magnitude(ctype: Optional[str], desc: dict) -> Optional[float]:
    if ctype == "mw":
        return max(0.0, desc["mw"] - 500)
    if ctype == "logp":
        return max(0.0, desc["logp"] - 5)
    if ctype == "hbd":
        return max(0.0, float(desc["hbd"] - 5))
    if ctype == "hba":
        return max(0.0, float(desc["hba"] - 10))
    if ctype == "rot":
        return max(0.0, float(desc["rot"] - 10))
    if ctype == "psa":
        return max(0.0, desc["psa"] - 140)
    if ctype in ("pains", "toxicity"):
        return 1.0
    return None


def _passed() -> dict:
    return {
        "passed": True,
        "failing_lemma": None,
        "failing_hypothesis": None,
        "constraint_type": None,
        "violation_magnitude": None,
    }
