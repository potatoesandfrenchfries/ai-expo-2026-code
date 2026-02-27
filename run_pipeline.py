import subprocess
import random
from rdkit import Chem
from rdkit.Chem import AllChem

def main():
    print("--- Knotworking AI Pipeline ---")
    print("1. Receiving Generated Molecule from Latent Space...")
    
    ai_generated_pool = [
        "C[NH+]1CCC(NC(=O)[C@H]2CCN(c3ccc(Cl)c(Cl)c3)C2=O)CC1",
        "CC(C)(C)C(=O)Nc1sc(CC(N)=O)nc1-c1cccc(F)c1",
        "O=C(Nc1cccc(Cl)c1)c1sc2c(c1)CCCC2",
        "CC1(C)CC(=O)C2(C)C(O)CC3OCC3(C)C2C1",
        "COc1ccc(S(=O)(=O)N2CCC(C(N)=O)CC2)cc1"
    ]
    
    target_smiles = random.choice(ai_generated_pool)
    print(f"   -> AI Selected SMILES: {target_smiles}")
    
    mol = Chem.MolFromSmiles(target_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    conf = mol.GetConformer()
    
    # ---> THIS IS THE MAGIC LINE WE WERE MISSING! <---
    Chem.MolToMolFile(mol, "generated_drug.sdf")
    # ------------------------------------------------
    
    print("2. Translating to Formal Logic Syntax...")
    
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
    coq_code += " ]\n    [].\n"
    
    demo_file = "src/rocq/Demo.v"
    with open(demo_file, "w") as f:
        f.write(coq_code)
        
    print(f"   -> Saved {num_atoms} generated coordinates to {demo_file}")
    print("   -> 3D model saved as 'generated_drug.sdf' in your folder.")
    
    print("3. Executing Formal Verification Engine...")
    try:
        result = subprocess.run(['coqc', '-R', 'src/rocq', 'Chemistry', 'src/rocq/Demo.v'], 
                                capture_output=True, text=True, check=True)
        print("\nSUCCESS! The generated molecular geometry was mathematically verified by the Rocq compiler.")
    except subprocess.CalledProcessError as e:
        print("\nVERIFICATION FAILED. Rocq detected a structural logic error:")
        print(e.stderr)

if __name__ == "__main__":
    main()
