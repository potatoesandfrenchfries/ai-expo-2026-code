import os

from rdkit import Chem
from rdkit.Chem import AllChem

from src.llm.pipeline_layer9 import Layer9FormalValidator
from src.llm.feedback_controller import FeedbackController


def main():
    print("--- Knotworking AI Pipeline ---")
    print("1. Running Generative Model (Bayesian VAE)...")

    target_smiles = "C[NH+]1CCC(NC(=O)[C@H]2CCN(c3ccc(Cl)c(Cl)c3)C2=O)CC1"
    mol = Chem.MolFromSmiles(target_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    print("2. Layer 9: Formal Validation via LLM-generated Rocq proof...")
    validator = Layer9FormalValidator()
    result = validator.validate_molecule(target_smiles)
    validator.print_report(result)

    # Persist the proof so it can be audited / used in Layer 22 IND drafting
    if result.proof:
        os.makedirs("data", exist_ok=True)
        proof_path = "data/last_proof.v"
        with open(proof_path, "w") as f:
            f.write(result.proof)
        print(f"\n  Rocq proof written → {proof_path}")

    # Layer 9 → Layer 7 feedback
    if not result.valid:
        print("\n3. Extracting generation constraints from failures...")
        controller = FeedbackController()
        constraints = controller.extract_constraints_from_failures([result])
        print(f"  Updated constraints: max_mw={constraints.max_mw}, "
              f"max_logP={constraints.max_logP}, "
              f"max_rot_bonds={constraints.max_rot_bonds}")
        print("  (Constraints will be applied to next VAE generation batch.)")


if __name__ == "__main__":
    main()
