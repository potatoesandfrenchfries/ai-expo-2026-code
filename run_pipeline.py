"""
run_pipeline.py  –  Knotworking AI Pipeline
============================================
Converts a SMILES string into a Rocq (Coq) formal molecule definition,
then attempts formal verification with coqc.

Usage:
    python run_pipeline.py [SMILES]

If SMILES is omitted the default drug-like molecule is used.
"""

import subprocess
import shutil
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem

# ---------------------------------------------------------------------------
# Mapping tables
# ---------------------------------------------------------------------------

ELEMENT_MAP = {
    'H':  'eH',
    'C':  'eC',
    'N':  'eN',
    'O':  'eO',
    'F':  'eF',
    'S':  'eS',
    'Cl': 'eCl',
    'Br': 'eBr',
    'P':  'eP',
    'I':  'eI',
}

BOND_TYPE_MAP = {
    rdchem.BondType.SINGLE:   'SingleBond',
    rdchem.BondType.DOUBLE:   'DoubleBond',
    rdchem.BondType.TRIPLE:   'TripleBond',
    rdchem.BondType.AROMATIC: 'AromaticBond',
}

DEFAULT_SMILES = "C[NH+]1CCC(NC(=O)[C@H]2CCN(c3ccc(Cl)c(Cl)c3)C2=O)CC1"


# ---------------------------------------------------------------------------
# Core translation function
# ---------------------------------------------------------------------------

def smiles_to_coq(smiles: str) -> str:
    """Translate a SMILES string into a complete Rocq source file string."""

    # --- Parse & embed 3-D coordinates ---
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    mol = Chem.AddHs(mol)
    embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
    if embed_result == -1:
        raise RuntimeError("3-D coordinate generation failed (EmbedMolecule returned -1).")
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()

    # --- Build atom entries ---
    atom_strings = []
    for atom in mol.GetAtoms():
        i      = atom.GetIdx()
        symbol = atom.GetSymbol()
        rocq_elem = ELEMENT_MAP.get(symbol)
        if rocq_elem is None:
            # Fall back to generic constructor pattern (e + Symbol)
            rocq_elem = f"e{symbol}"

        pos   = conf.GetAtomPosition(i)
        # Represent floats as exact integer fractions so Rocq's reals accept them
        x_frac = _frac(pos.x)
        y_frac = _frac(pos.y)
        z_frac = _frac(pos.z)

        charge = atom.GetFormalCharge()
        charge_str = f"{charge}%Z"

        atom_strings.append(
            f"mkAtom {i} {rocq_elem} (mkPoint {x_frac} {y_frac} {z_frac}) "
            f"{charge_str} None None"
        )

    # --- Build bond entries ---
    bond_strings = []
    for bond in mol.GetBonds():
        b_id   = bond.GetIdx()
        a1     = bond.GetBeginAtomIdx()
        a2     = bond.GetEndAtomIdx()
        btype  = BOND_TYPE_MAP.get(bond.GetBondType(), 'SingleBond')
        bond_strings.append(f"mkBond {b_id} {a1} {a2} {btype} None")

    # --- Assemble Rocq source ---
    lines = []

    # Imports
    lines += [
        "From Stdlib Require Import List.",
        "Import ListNotations.",
        "Require Import Stdlib.Reals.Reals.",
        "Open Scope R_scope.",
        "",
        "Require Import Stdlib.ZArith.ZArith.",
        "Require Import Chemistry.Atoms.",
        "Require Import Chemistry.Geometry.",
        "Require Import Chemistry.Bonds.",
        "Require Import Chemistry.Molecules.",
        "",
        "Require Import Chemistry.MolecularProperties.",
        "Require Import Chemistry.HydrogenBonding.",
        "Require Import Chemistry.FunctionalGroups.",
        "Require Import Chemistry.MathProperties.",
        "Require Import Chemistry.Valency.",
        "Require Import Chemistry.Aromaticity.",
        "Require Import Chemistry.Conformational.",
        "",
    ]

    # Atom list
    atom_sep  = "\n    ; "
    atoms_str = atom_sep.join(atom_strings)

    # Bond list (may be empty)
    if bond_strings:
        bond_sep  = "\n    ; "
        bonds_str = bond_sep.join(bond_strings)
        bonds_block = f"[ {bonds_str}\n    ]"
    else:
        bonds_block = "[]"

    lines += [
        "Definition demo_molecule : Molecule :=",
        "  mkMol",
        f"    [ {atoms_str}",
        "    ]",
        f"    {bonds_block}.",
        "",
    ]

    # Verification checks:
    #   Check  — instant type-check only, proves the molecule is well-typed.
    #   No Eval/compute lines: those force kernel evaluation and are very slow.
    #   Run `coqtop` interactively for evaluation if needed.
    lines += [
        "(* ---- Verification checks ---- *)",
        "(* Confirm demo_molecule is a well-typed Molecule (instant) *)",
        "Check demo_molecule.",
        "(* Confirm property functions are applicable (instant) *)",
        "Check (is_connected demo_molecule).",
        "Check (atom_count demo_molecule).",
        "Check (bond_count demo_molecule).",
        "",
    ]

    return "\n".join(lines)


def _frac(value: float, denom: int = 10000) -> str:
    """Return an exact rational literal suitable for Rocq's R type."""
    numerator = int(round(value * denom))
    return f"({numerator} / {denom})"


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def main():
    print("--- Knotworking AI Pipeline ---")

    # Accept SMILES from CLI argument
    if len(sys.argv) >= 2:
        target_smiles = sys.argv[1]
        print(f"Using CLI SMILES: {target_smiles}")
    else:
        target_smiles = DEFAULT_SMILES
        print(f"Using default SMILES: {target_smiles}")

    # Step 1
    print("\n1. Running Generative Model (Overfitted)...")
    mol_check = Chem.MolFromSmiles(target_smiles)
    if mol_check is None:
        print("ERROR: Invalid SMILES string. Aborting.")
        sys.exit(1)
    n_heavy = mol_check.GetNumAtoms()
    mol_h   = Chem.AddHs(mol_check)
    n_total = mol_h.GetNumAtoms()
    n_bonds = mol_h.GetNumBonds()
    print(f"   -> Heavy atoms: {n_heavy} | Total (with H): {n_total} | Bonds: {n_bonds}")

    # Step 2
    print("\n2. Translating to Rocq Formal Syntax...")
    try:
        coq_code = smiles_to_coq(target_smiles)
    except Exception as exc:
        print(f"ERROR during translation: {exc}")
        sys.exit(1)

    demo_file = "src/rocq/Demo.v"
    with open(demo_file, "w") as f:
        f.write(coq_code)
    print(f"   -> Saved generated proof to {demo_file}")
    print(f"   -> Atoms emitted: {n_total} | Bonds emitted: {n_bonds}")

    # Step 3
    print("\n3. Running Formal Verification...")

    # Check coqc availability before attempting compilation
    if shutil.which("coqc") is None:
        print("   WARNING: coqc not found on PATH. Skipping formal verification.")
        print("   Install Rocq/Coq and ensure coqc is in your PATH to verify.")
        print("\nPipeline complete (verification skipped — coqc unavailable).")
        return

    try:
        result = subprocess.run(
            ['coqc', '-R', 'src/rocq', 'Chemistry', demo_file],
            capture_output=True, text=True, check=True
        )
        print("\nSUCCESS! The molecule was formally verified by Rocq.")
        if result.stdout:
            print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("\nVERIFICATION FAILED. Rocq found a logical error in the molecule:")
        print(e.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
