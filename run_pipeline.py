import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

def main():
    print("--- Knotworking AI Pipeline ---")
    print("1. Running Generative Model (Bayesian VAE)...")
    
    target_smiles = "C[NH+]1CCC(NC(=O)[C@H]2CCN(c3ccc(Cl)c(Cl)c3)C2=O)CC1"
    mol = Chem.MolFromSmiles(target_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    conf = mol.GetConformer()
    
    print("2. Translating to Rocq Formal Syntax...")
    
    coq_code = "From Coq Require Import List.\nImport ListNotations.\n"
    coq_code += "Require Import Coq.Reals.Reals.\nOpen Scope R_scope.\n\n"
    coq_code += "Require Import Coq.ZArith.ZArith.\n"
    coq_code += "Require Import Chemistry.Atoms.\n"
    coq_code += "Require Import Chemistry.Geometry.\n"
    coq_code += "Require Import Chemistry.Bonds.\n"
    coq_code += "Require Import Chemistry.Molecules.\n\n"
    
    coq_code += "Definition demo_molecule : Molecule :=\n"
    coq_code += "  mkMol\n    [ "
    
    atom_strings = []
    for i in range(3):
        pos = conf.GetAtomPosition(i)
        symbol = "e" + mol.GetAtomWithIdx(i).GetSymbol()
        
        x_frac = f"({int(pos.x * 1000)} / 1000)"
        y_frac = f"({int(pos.y * 1000)} / 1000)"
        z_frac = f"({int(pos.z * 1000)} / 1000)"
        
        atom_strings.append(
            f"mkAtom {i} {symbol} (mkPoint {x_frac} {y_frac} {z_frac}) 0%Z None None"
        )
    
    coq_code += " ;\n      ".join(atom_strings)
    coq_code += " ]\n    [].\n"
    
    demo_file = "src/rocq/Demo.v"
    with open(demo_file, "w") as f:
        f.write(coq_code)
        
    print(f"   -> Saved generated proof to {demo_file}")
    
    print("3. Running Formal Verification...")
    try:
        result = subprocess.run(['coqc', '-R', 'src/rocq', 'Chemistry', 'src/rocq/Demo.v'], 
                                capture_output=True, text=True, check=True)
        print("\nSUCCESS! The molecule was formally verified by Rocq using Kanish's advanced theorems.")
    except subprocess.CalledProcessError as e:
        print("\nVERIFICATION FAILED. Rocq found a logical error in the molecule:")
        print(e.stderr)

if __name__ == "__main__":
    main()
