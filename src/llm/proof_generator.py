"""LLM-based Rocq proof generator using the Claude API (claude-haiku-4-5).

Converts a SMILES string + RDKit molecular descriptors into a formally
verifiable Rocq proof of drug-likeness (Lipinski RO5 by default).
Includes a retry loop that feeds coqc error output back to the model.
"""

from dataclasses import dataclass
from typing import Optional

import anthropic

from .rocq_verifier import RocqVerifier, VerificationResult

_SYSTEM_PROMPT = """\
You are a Rocq/Coq formal proof assistant specialising in computational chemistry.
You generate compilable Rocq source files that prove drug-likeness properties
of small molecules using the Chemistry library.

Key types available after imports:
  MolDescriptors  – record with fields:
      md_mol_weight  md_logP  md_hbd  md_hba  md_rot_bonds
      md_psa  md_molar_refract  md_logS  md_rings
      md_arom_rings  md_chiral_centers  md_atom_count
  lipinski_ro5 : MolDescriptors -> Prop   (MW≤500, logP≤5, HBD≤5, HBA≤10)
  veber_rules  : MolDescriptors -> Prop   (rot_bonds≤10, PSA≤140)

Standard proof tactics:
  - Real inequalities  → lra
  - Nat inequalities   → lia
  - Conjunctions       → split. { lra. } split. { lra. } ...

Rules:
1. Always start the file with the required imports block shown below.
2. Define a MolDescriptors record named mol_desc using mkMolDesc.
3. State and prove (or disprove) the requested lemma.
4. Respond ONLY with valid Rocq code — no prose, no markdown fences.

Required imports block:
From Coq Require Import List.
Import ListNotations.
Require Import Coq.Reals.Reals.
Open Scope R_scope.
Require Import Coq.ZArith.ZArith.
Require Import Chemistry.DrugLikeness.
Require Import Coq.micromega.Lra.
Require Import Coq.micromega.Lia.
"""

# Two verified few-shot examples that ground the model in correct syntax
_FEW_SHOT: list[dict] = [
    {
        "role": "user",
        "content": (
            "SMILES: CC1=CC=CC=C1\nMW: 92.14\nLogP: 2.73\nHBD: 0\nHBA: 0\n"
            "RotBonds: 0\nPSA: 0.0\nMolRefract: 31.27\nLogS: -2.0\n"
            "Rings: 1\nAromRings: 1\nChiralCenters: 0\nHeavyAtoms: 7\n"
            "Task: Prove Lipinski RO5"
        ),
    },
    {
        "role": "assistant",
        "content": (
            "From Coq Require Import List.\n"
            "Import ListNotations.\n"
            "Require Import Coq.Reals.Reals.\n"
            "Open Scope R_scope.\n"
            "Require Import Coq.ZArith.ZArith.\n"
            "Require Import Chemistry.DrugLikeness.\n"
            "Require Import Coq.micromega.Lra.\n"
            "Require Import Coq.micromega.Lia.\n\n"
            "Definition mol_desc : MolDescriptors :=\n"
            "  mkMolDesc 92.14 2.73 0 0 0 0.0 31.27 (-2.0) 1 1 0 7.\n\n"
            "Lemma mol_lipinski : lipinski_ro5 mol_desc.\n"
            "Proof.\n"
            "  unfold lipinski_ro5, lipinski_mw, lipinski_logP,\n"
            "         lipinski_hbd, lipinski_hba, mol_desc.\n"
            "  simpl. split. { lra. } split. { lra. } split. { lia. } { lia. }\n"
            "Qed.\n"
        ),
    },
    {
        "role": "user",
        "content": (
            "SMILES: c1ccc(cc1)C(=O)O\nMW: 122.12\nLogP: 1.87\nHBD: 1\nHBA: 2\n"
            "RotBonds: 1\nPSA: 37.3\nMolRefract: 32.89\nLogS: -1.56\n"
            "Rings: 1\nAromRings: 1\nChiralCenters: 0\nHeavyAtoms: 9\n"
            "Task: Prove Lipinski RO5"
        ),
    },
    {
        "role": "assistant",
        "content": (
            "From Coq Require Import List.\n"
            "Import ListNotations.\n"
            "Require Import Coq.Reals.Reals.\n"
            "Open Scope R_scope.\n"
            "Require Import Coq.ZArith.ZArith.\n"
            "Require Import Chemistry.DrugLikeness.\n"
            "Require Import Coq.micromega.Lra.\n"
            "Require Import Coq.micromega.Lia.\n\n"
            "Definition mol_desc : MolDescriptors :=\n"
            "  mkMolDesc 122.12 1.87 1 2 1 37.3 32.89 (-1.56) 1 1 0 9.\n\n"
            "Lemma mol_lipinski : lipinski_ro5 mol_desc.\n"
            "Proof.\n"
            "  unfold lipinski_ro5, lipinski_mw, lipinski_logP,\n"
            "         lipinski_hbd, lipinski_hba, mol_desc.\n"
            "  simpl. split. { lra. } split. { lra. } split. { lia. } { lia. }\n"
            "Qed.\n"
        ),
    },
]


@dataclass
class ProofResult:
    success: bool
    proof: str
    attempts: int
    error_type: Optional[str]


class RocqProofGenerator:
    """Generates and iteratively repairs Rocq proofs via the Claude API."""

    def __init__(self, verifier: Optional[RocqVerifier] = None):
        self.client = anthropic.Anthropic()
        self.verifier = verifier or RocqVerifier()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def generate(
        self,
        smiles: str,
        descriptors: dict,
        task: str = "Prove Lipinski RO5",
    ) -> str:
        """Single-shot generation with no verification."""
        messages = _FEW_SHOT + [
            {"role": "user", "content": self._build_prompt(smiles, descriptors, task)}
        ]
        response = self.client.messages.create(
            model="claude-haiku-4-5-20251001",
            max_tokens=1024,
            temperature=0.1,
            system=_SYSTEM_PROMPT,
            messages=messages,
        )
        return response.content[0].text.strip()

    def generate_with_retry(
        self,
        smiles: str,
        descriptors: dict,
        task: str = "Prove Lipinski RO5",
        max_attempts: int = 3,
    ) -> ProofResult:
        """Generate, verify with coqc, and repair on failure (up to max_attempts)."""
        messages = list(_FEW_SHOT) + [
            {"role": "user", "content": self._build_prompt(smiles, descriptors, task)}
        ]
        proof_text = ""
        verification: Optional[VerificationResult] = None

        for attempt in range(1, max_attempts + 1):
            response = self.client.messages.create(
                model="claude-haiku-4-5-20251001",
                max_tokens=1024,
                temperature=0.1,
                system=_SYSTEM_PROMPT,
                messages=messages,
            )
            proof_text = response.content[0].text.strip()
            verification = self.verifier.verify(proof_text)

            if verification.success:
                return ProofResult(
                    success=True,
                    proof=proof_text,
                    attempts=attempt,
                    error_type=None,
                )

            # Append model output + correction request for next round
            messages.append({"role": "assistant", "content": proof_text})
            error_summary = "\n".join(verification.error_lines[:8])
            messages.append(
                {
                    "role": "user",
                    "content": (
                        f"That proof failed to compile "
                        f"(error type: {verification.error_type}).\n"
                        f"coqc output:\n{error_summary}\n\n"
                        "Please correct the proof. Return only valid Rocq code."
                    ),
                }
            )

        return ProofResult(
            success=False,
            proof=proof_text,
            attempts=max_attempts,
            error_type=verification.error_type if verification else "UNKNOWN",
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _build_prompt(smiles: str, descriptors: dict, task: str) -> str:
        lines = [f"SMILES: {smiles}"]
        for key, val in descriptors.items():
            lines.append(f"{key}: {val}")
        lines.append(f"Task: {task}")
        return "\n".join(lines)
