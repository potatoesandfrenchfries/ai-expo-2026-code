"""Training data generator for the Rocq proof LLM.

Reads molecules_clean.csv (produced by data_loader.py), extracts RDKit
descriptors, templates a compilable Rocq Lipinski proof for each molecule,
verifies it with coqc, and writes verified (input, output) pairs to
data/rocq_proof_corpus.jsonl for fine-tuning or few-shot retrieval.

Usage:
    python -m src.llm.training_data_generator
"""

import json
import os

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from .rocq_verifier import RocqVerifier

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "data")
OUTPUT_PATH = os.path.join(DATA_DIR, "rocq_proof_corpus.jsonl")

_IMPORTS = (
    "From Coq Require Import List.\n"
    "Import ListNotations.\n"
    "Require Import Coq.Reals.Reals.\n"
    "Open Scope R_scope.\n"
    "Require Import Coq.ZArith.ZArith.\n"
    "Require Import Chemistry.DrugLikeness.\n"
    "Require Import Coq.micromega.Lra.\n"
    "Require Import Coq.micromega.Lia.\n\n"
)


# ---------------------------------------------------------------------------
# Descriptor extraction
# ---------------------------------------------------------------------------

def extract_descriptors(mol) -> dict:
    """Return a flat dict of RDKit descriptors keyed by display name."""
    logP = Descriptors.MolLogP(mol)
    return {
        "MW":            round(Descriptors.ExactMolWt(mol), 4),
        "LogP":          round(logP, 4),
        "HBD":           rdMolDescriptors.CalcNumHBD(mol),
        "HBA":           rdMolDescriptors.CalcNumHBA(mol),
        "RotBonds":      rdMolDescriptors.CalcNumRotatableBonds(mol),
        "PSA":           round(Descriptors.TPSA(mol), 4),
        "MolRefract":    round(Descriptors.MolMR(mol), 4),
        # Moriguchi logS estimate
        "LogS":          round(-logP * 0.54 - 0.54, 4),
        "Rings":         rdMolDescriptors.CalcNumRings(mol),
        "AromRings":     rdMolDescriptors.CalcNumAromaticRings(mol),
        "ChiralCenters": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
        "HeavyAtoms":    mol.GetNumHeavyAtoms(),
    }


# ---------------------------------------------------------------------------
# Proof templating
# ---------------------------------------------------------------------------

def _passes_lipinski(d: dict) -> bool:
    return (
        d["MW"] <= 500
        and d["LogP"] <= 5
        and d["HBD"] <= 5
        and d["HBA"] <= 10
    )


def build_lipinski_proof(smiles: str, d: dict) -> str:
    """Template a Rocq file that proves (or refutes) lipinski_ro5 for d."""
    logS_str = f"({d['LogS']})" if d["LogS"] < 0 else str(d["LogS"])

    proof = _IMPORTS
    proof += (
        f"Definition mol_desc : MolDescriptors :=\n"
        f"  mkMolDesc {d['MW']} {d['LogP']} {d['HBD']} {d['HBA']} "
        f"{d['RotBonds']} {d['PSA']} {d['MolRefract']} {logS_str} "
        f"{d['Rings']} {d['AromRings']} {d['ChiralCenters']} {d['HeavyAtoms']}.\n\n"
    )

    if _passes_lipinski(d):
        proof += (
            "Lemma mol_lipinski : lipinski_ro5 mol_desc.\n"
            "Proof.\n"
            "  unfold lipinski_ro5, lipinski_mw, lipinski_logP,\n"
            "         lipinski_hbd, lipinski_hba, mol_desc.\n"
            "  simpl. split. { lra. } split. { lra. } split. { lia. } { lia. }\n"
            "Qed.\n"
        )
    else:
        proof += (
            "Lemma mol_not_lipinski : ~ lipinski_ro5 mol_desc.\n"
            "Proof.\n"
            "  unfold lipinski_ro5, lipinski_mw, lipinski_logP,\n"
            "         lipinski_hbd, lipinski_hba, mol_desc.\n"
            "  simpl. intro H.\n"
            "  destruct H as [Hmw [HlogP [Hhbd Hhba]]]. lra.\n"
            "Qed.\n"
        )
    return proof


def build_input_prompt(smiles: str, d: dict, task: str) -> str:
    lines = [f"SMILES: {smiles}"]
    for k, v in d.items():
        lines.append(f"{k}: {v}")
    lines.append(f"Task: {task}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main corpus generation
# ---------------------------------------------------------------------------

def generate_corpus(max_molecules: int = 500) -> None:
    """
    Generate the Rocq proof corpus and write it to data/rocq_proof_corpus.jsonl.

    Each entry:
        smiles       – canonical SMILES
        input        – prompt text (descriptors + task)
        output       – templated Rocq proof source
        verified     – True if coqc compiled successfully
        error_type   – error taxonomy string or None
    """
    csv_path = os.path.join(DATA_DIR, "molecules_clean.csv")
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"{csv_path} not found — run src/data/data_loader.py first."
        )

    df = pd.read_csv(csv_path).head(max_molecules)
    verifier = RocqVerifier()
    corpus: list[dict] = []
    verified_count = 0

    for _, row in df.iterrows():
        smiles = row["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        d = extract_descriptors(mol)
        task = "Prove Lipinski RO5"
        proof_text = build_lipinski_proof(smiles, d)
        result = verifier.verify(proof_text)

        if result.success:
            verified_count += 1

        corpus.append(
            {
                "smiles":     smiles,
                "input":      build_input_prompt(smiles, d, task),
                "output":     proof_text,
                "verified":   result.success,
                "error_type": result.error_type,
            }
        )

    os.makedirs(DATA_DIR, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        for entry in corpus:
            f.write(json.dumps(entry) + "\n")

    total = len(corpus)
    print(f"Corpus complete: {total} molecules, {verified_count} verified proofs "
          f"({100 * verified_count / max(total, 1):.1f}%)")
    print(f"Saved → {OUTPUT_PATH}")


if __name__ == "__main__":
    generate_corpus()
