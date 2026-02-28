"""
Knotworking AI — Verified Drug Discovery Pipeline

    1. Generate/select a candidate molecule (Layer 7)
    2. Generate 3D coordinates with RDKit
    3. Verify molecular structure with Rocq (Layer 9a)
    4. Encode molecule as fingerprint and predict properties with BayesianGraphVAE (Layer 8)
    5. Verify predictions and drug-likeness with Rocq (Layer 9b)
"""

import subprocess
import random
import os
import torch
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs

from model import BayesianGraphVAE, mc_predict


OUTPUT_DIR = "outputs"
ROCQ_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "rocq")
MODEL_PATH = "models/model.pt"


def ensure_output_dirs():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(ROCQ_OUTPUT_DIR, exist_ok=True)


# ──────────────────────────────────────────────
# MOLECULE → FINGERPRINT
# ──────────────────────────────────────────────

def molecule_to_fingerprint(mol, n_bits=2048):
    """Convert an RDKit molecule to a 2048-bit Morgan fingerprint tensor."""
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
    arr = np.zeros(n_bits, dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return torch.tensor(arr).unsqueeze(0)


def rdkit_bond_to_rocq(rdkit_bond_type):
    if rdkit_bond_type == Chem.rdchem.BondType.SINGLE:
        return "SingleBond"
    if rdkit_bond_type == Chem.rdchem.BondType.DOUBLE:
        return "DoubleBond"
    if rdkit_bond_type == Chem.rdchem.BondType.TRIPLE:
        return "TripleBond"
    if rdkit_bond_type == Chem.rdchem.BondType.AROMATIC:
        return "AromaticBond"
    return "SingleBond"


# ──────────────────────────────────────────────
# ROCQ CODE GENERATION & COMPILATION
# ──────────────────────────────────────────────

def generate_structure_rocq(mol, conf):
    """Generate Rocq code to verify molecular structure (valid atoms, geometry)."""

    coq_code = "From Stdlib Require Import List.\nImport ListNotations.\n"
    coq_code += "Require Import Stdlib.Reals.Reals.\nOpen Scope R_scope.\n\n"
    coq_code += "Require Import Stdlib.ZArith.ZArith.\n"
    coq_code += "Require Import Stdlib.Arith.Arith.\n"
    coq_code += "Require Import Stdlib.micromega.Lia.\n"
    coq_code += "Require Import Chemistry.Atoms.\n"
    coq_code += "Require Import Chemistry.Geometry.\n"
    coq_code += "Require Import Chemistry.Bonds.\n"
    coq_code += "Require Import Chemistry.Molecules.\n\n"

    coq_code += "Definition demo_molecule : Molecule :=\n"
    coq_code += "  mkMol\n    [ "

    atom_strings = []
    num_atoms = mol.GetNumAtoms()
    for i in range(num_atoms):
        pos = conf.GetAtomPosition(i)
        symbol = "e" + mol.GetAtomWithIdx(i).GetSymbol()

        x_frac = f"({int(pos.x * 1000)} / 1000)"
        y_frac = f"({int(pos.y * 1000)} / 1000)"
        z_frac = f"({int(pos.z * 1000)} / 1000)"

        atom_strings.append(
            f"mkAtom {i} {symbol} (mkPoint {x_frac} {y_frac} {z_frac}) 0%Z None None"
        )

    coq_code += " ;\n      ".join(atom_strings)
    coq_code += " ]\n    [ "

    bond_strings = []
    for i, bond in enumerate(mol.GetBonds()):
        bond_strings.append(
            f"mkBond {i} {bond.GetBeginAtomIdx()} {bond.GetEndAtomIdx()} {rdkit_bond_to_rocq(bond.GetBondType())} None"
        )

    coq_code += " ;\n      ".join(bond_strings)
    coq_code += " ].\n\n"

    coq_code += "Theorem generated_has_atoms : length (mol_atoms demo_molecule) > 0.\n"
    coq_code += "Proof. compute; lia. Qed.\n\n"
    coq_code += "Theorem generated_has_bonds : length (mol_bonds demo_molecule) > 0.\n"
    coq_code += "Proof. compute; lia. Qed.\n"

    return coq_code


def generate_decision_rocq(mean_var, var_spread, mol):
    """Generate Rocq code to verify BNN predictions meet drug-likeness criteria."""

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    coq_code = "Require Import Stdlib.Reals.Reals.\n"
    coq_code += "Require Import Stdlib.micromega.Lra.\n"
    coq_code += "Open Scope R_scope.\n\n"

    coq_code += "(* === BNN Uncertainty Values === *)\n"
    coq_code += f"Definition mean_variance : R := {mean_var:.6f}.\n"
    coq_code += f"Definition variance_spread : R := {var_spread:.6f}.\n\n"

    coq_code += "(* === Lipinski Molecular Descriptors (computed by RDKit) === *)\n"
    coq_code += f"Definition mol_weight : R := {mw:.4f}.\n"
    coq_code += f"Definition mol_logP : R := {logp:.4f}.\n"
    coq_code += f"Definition mol_hbd : R := {float(hbd):.1f}.\n"
    coq_code += f"Definition mol_hba : R := {float(hba):.1f}.\n\n"

    coq_code += "(* === Formal Verification Theorems === *)\n\n"

    coq_code += "(* Model reconstruction uncertainty must be below threshold *)\n"
    coq_code += "Theorem confidence_check : mean_variance < 0.5.\n"
    coq_code += "Proof. unfold mean_variance. lra. Qed.\n\n"

    coq_code += "(* Uncertainty spread must also remain bounded *)\n"
    coq_code += "Theorem spread_check : variance_spread < 0.5.\n"
    coq_code += "Proof. unfold variance_spread. lra. Qed.\n\n"

    coq_code += "(* Lipinski Rule of Five: molecular weight <= 500 *)\n"
    coq_code += "Theorem lipinski_mw : mol_weight <= 500.0.\n"
    coq_code += "Proof. unfold mol_weight. lra. Qed.\n\n"

    coq_code += "(* Lipinski Rule of Five: LogP <= 5 *)\n"
    coq_code += "Theorem lipinski_logp : mol_logP <= 5.0.\n"
    coq_code += "Proof. unfold mol_logP. lra. Qed.\n\n"

    coq_code += "(* Lipinski Rule of Five: hydrogen bond donors <= 5 *)\n"
    coq_code += "Theorem lipinski_hbd : mol_hbd <= 5.0.\n"
    coq_code += "Proof. unfold mol_hbd. lra. Qed.\n\n"

    coq_code += "(* Lipinski Rule of Five: hydrogen bond acceptors <= 10 *)\n"
    coq_code += "Theorem lipinski_hba : mol_hba <= 10.0.\n"
    coq_code += "Proof. unfold mol_hba. lra. Qed.\n\n"

    return coq_code


def run_rocq(filepath):
    """Compile a .v file with coqc. Returns (success, error_message)."""
    try:
        result = subprocess.run(
            ['coqc', '-R', 'src/rocq', 'Chemistry', filepath],
            capture_output=True, text=True, check=True
        )
        return True, ""
    except subprocess.CalledProcessError as e:
        return False, e.stderr


# ──────────────────────────────────────────────
# MAIN PIPELINE
# ──────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  KNOTWORKING AI — Verified Drug Discovery Pipeline")
    print("=" * 60)

    ensure_output_dirs()


    ### 1. Simulate receiving a molecule from the latent space of a generative model trained on drug-like molecules
    print("\n[Layer 7] Generating candidate molecule...")

    # Hardcoded pool of molecules - ideally this would come from the trained BayesianGraphVAE's decoder
    ai_generated_pool = [
        "C[NH+]1CCC(NC(=O)[C@H]2CCN(c3ccc(Cl)c(Cl)c3)C2=O)CC1",
        "CC(C)(C)C(=O)Nc1sc(CC(N)=O)nc1-c1cccc(F)c1",
        "O=C(Nc1cccc(Cl)c1)c1sc2c(c1)CCCC2",
        "CC1(C)CC(=O)C2(C)C(O)CC3OCC3(C)C2C1",
        "COc1ccc(S(=O)(=O)N2CCC(C(N)=O)CC2)cc1"
    ]

    target_smiles = random.choice(ai_generated_pool)
    print(f"   SMILES: {target_smiles}")


    ### 2. Generate 3D coordinates for the selected molecule using RDKit
    mol = Chem.MolFromSmiles(target_smiles)
    mol_with_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_with_h, randomSeed=42)
    conf = mol_with_h.GetConformer()
    sdf_path = os.path.join(OUTPUT_DIR, "generated_drug.sdf")
    Chem.MolToMolFile(mol_with_h, sdf_path)
    print(f"   Atoms: {mol_with_h.GetNumAtoms()} | Bonds: {mol_with_h.GetNumBonds()}")
    print(f"   3D model saved as '{sdf_path}'")


    ### 3. Verify molecular structure with Rocq (Layer 9a)
    print("\nVerifying molecular structure...")

    structure_code = generate_structure_rocq(mol_with_h, conf)
    structure_file = os.path.join(ROCQ_OUTPUT_DIR, "GeneratedDemo.v")
    with open(structure_file, "w") as f:
        f.write(structure_code)

    success, error = run_rocq(structure_file)
    if success:
        print("   ✓ Rocq compilation succeeded for generated structure module")
    else:
        print("   ✗ Structure verification FAILED:")
        print(f"     {error}")
        print("   Pipeline halted — molecule is structurally invalid.")
        return


    ### 4. Encode molecule as fingerprint and predict properties with BayesianGraphVAE (Layer 8)
    print("\n[Layer 8] Predicting molecular properties...")

    # Convert molecule to 2048-bit Morgan fingerprint (the format the VAE was trained on)
    fingerprint = molecule_to_fingerprint(mol)

    # Load the trained model
    if not os.path.exists(MODEL_PATH):
        print(f"   ✗ Missing trained model checkpoint: '{MODEL_PATH}'")
        print("   Run 'python model.py' first to train and save the model.")
        return

    model = BayesianGraphVAE(input_dim=2048)
    model.load_state_dict(torch.load(MODEL_PATH, map_location="cpu", weights_only=True))

    # Monte Carlo uncertainty estimation — 50 forward passes through the VAE
    # Each pass samples different latent vectors due to reparameterization
    # Variance across passes = model's uncertainty about this molecule
    mean_pred, variance = mc_predict(model, fingerprint, samples=50)

    mean_var = variance.mean().item()       # average uncertainty across all fingerprint bits
    var_spread = variance.std().item()      # how uneven the uncertainty is across bits
    recon_error = ((mean_pred - fingerprint) ** 2).mean().item()

    print(f"   Reconstruction error:   {recon_error:.6f}")
    print(f"   Mean variance:          {mean_var:.6f}")
    print(f"   Variance spread:        {var_spread:.6f}")


    ### 5. Verify predictions and drug-likeness with Rocq (Layer 9b)
    print("\n Verifying drug-likeness and prediction confidence...")

    decision_code = generate_decision_rocq(mean_var, var_spread, mol)
    decision_file = os.path.join(ROCQ_OUTPUT_DIR, "GeneratedDecision.v")
    with open(decision_file, "w") as f:
        f.write(decision_code)

    success, error = run_rocq(decision_file)
    if success:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        print("   ✓ Confidence check passed")
        print("   ✓ Lipinski Rule of Five satisfied")
        print(f"\n{'=' * 60}")
        print(f"  MOLECULE ACCEPTED")
        print(f"  SMILES:      {target_smiles}")
        print(f"  Uncertainty: {mean_var:.6f}")
        print(f"  MW: {mw:.1f} | LogP: {logp:.2f} | HBD: {hbd} | HBA: {hba}")
        print(f"{'=' * 60}")
    else:
        print("   ✗ Decision verification FAILED:")
        print(f"     {error}")
        print(f"\n  MOLECULE REJECTED — failed formal property verification")


if __name__ == "__main__":
    main()
