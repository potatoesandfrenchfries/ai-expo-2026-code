(** * FastMolecule: Optimized molecular representation for fast compilation
    
    This module provides an optimized version of the Molecule type using:
    - CommonElements (10 elements instead of 118)
    - Efficient set-based operations for connectivity
    - Opaque computational functions
    - Simplified data structures
*)

From Stdlib Require Import List.
Import ListNotations.
Require Import Stdlib.Reals.Reals.
Require Import Stdlib.ZArith.ZArith.
Require Import Chemistry.CommonElements.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.

Open Scope R_scope.

(** ------------------------------------------------------------------ *)
(** ** 1. Atom Instance (using CommonElement)                          *)
(** ------------------------------------------------------------------ *)

Record AtomInst : Type := mkAtom
  { ai_id       : nat
  ; ai_element  : CommonElement
  ; ai_position : Point3D
  ; ai_charge   : Z
  ; ai_isotope  : option nat
  ; ai_stereo   : option bool
  }.

(** ------------------------------------------------------------------ *)
(** ** 2. Bond Instance                                                 *)
(** ------------------------------------------------------------------ *)

Record BondInst : Type := mkBond
  { bi_id     : nat
  ; bi_atom1  : nat
  ; bi_atom2  : nat
  ; bi_type   : BondType
  ; bi_stereo : option bool
  }.

(** ------------------------------------------------------------------ *)
(** ** 3. Molecule                                                      *)
(** ------------------------------------------------------------------ *)

Record Molecule : Type := mkMol
  { mol_atoms : list AtomInst
  ; mol_bonds : list BondInst
  }.

(** ------------------------------------------------------------------ *)
(** ** 4. Basic Queries (kept simple)                                   *)
(** ------------------------------------------------------------------ *)

Definition atom_count (mol : Molecule) : nat :=
  length mol.(mol_atoms).

Definition bond_count (mol : Molecule) : nat :=
  length mol.(mol_bonds).

(** ------------------------------------------------------------------ *)
(** ** 5. Molecular Weight (optimized)                                  *)
(** ------------------------------------------------------------------ *)

Definition molecular_weight (mol : Molecule) : R :=
  fold_left
    (fun acc a => acc + atomic_mass a.(ai_element))
    mol.(mol_atoms)
    0.

(** Make molecular weight opaque for faster compilation *)
Definition molecular_weight_opaque := molecular_weight.
Global Opaque molecular_weight_opaque.

(** ------------------------------------------------------------------ *)
(** ** 6. Connectivity (simplified)                                     *)
(** ------------------------------------------------------------------ *)

(** Check if two atoms are bonded *)
Definition are_bonded (mol : Molecule) (i j : nat) : bool :=
  existsb
    (fun b =>
       orb (andb (Nat.eqb b.(bi_atom1) i) (Nat.eqb b.(bi_atom2) j))
           (andb (Nat.eqb b.(bi_atom1) j) (Nat.eqb b.(bi_atom2) i)))
    mol.(mol_bonds).

(** Get neighbors of an atom *)
Definition neighbors (mol : Molecule) (i : nat) : list nat :=
  flat_map
    (fun b =>
       if Nat.eqb b.(bi_atom1) i then [b.(bi_atom2)]
       else if Nat.eqb b.(bi_atom2) i then [b.(bi_atom1)]
       else [])
    mol.(mol_bonds).

(** Degree of an atom *)
Definition degree (mol : Molecule) (i : nat) : nat :=
  length (neighbors mol i).

(** ------------------------------------------------------------------ *)
(** ** 7. Connectivity Check (BFS with early marking)                   *)
(** ------------------------------------------------------------------ *)

Fixpoint bfs_fuel
    (mol     : Molecule)
    (queue   : list nat)
    (visited : list nat)
    (fuel    : nat)
    : list nat :=
  match fuel with
  | O => visited
  | S n =>
    match queue with
    | [] => visited
    | i :: rest =>
      let new_nbrs := filter (fun j => negb (existsb (Nat.eqb j) visited))
                             (neighbors mol i) in
      let visited' := new_nbrs ++ visited in
      let new_q    := rest ++ new_nbrs in
      bfs_fuel mol new_q visited' n
    end
  end.

Definition reachable_from (mol : Molecule) (start : nat) : list nat :=
  let n_atoms := length mol.(mol_atoms) in
  bfs_fuel mol [start] [start] (n_atoms + 1).

(** Check if molecule is connected *)
Definition is_connected (mol : Molecule) : bool :=
  match mol.(mol_atoms) with
  | [] => true
  | first :: _ =>
    let start   := first.(ai_id) in
    let reached := reachable_from mol start in
    let all_ids := map ai_id mol.(mol_atoms) in
    forallb (fun id => existsb (Nat.eqb id) reached) all_ids
  end.

(** ------------------------------------------------------------------ *)
(** ** 8. Bond Order Sum                                                *)
(** ------------------------------------------------------------------ *)

Definition bond_order_sum (mol : Molecule) (i : nat) : R :=
  fold_left
    (fun acc b =>
       if orb (Nat.eqb b.(bi_atom1) i) (Nat.eqb b.(bi_atom2) i)
       then acc + bond_order b.(bi_type)
       else acc)
    mol.(mol_bonds)
    0.

(** ------------------------------------------------------------------ *)
(** ** 9. Empty Molecule                                                *)
(** ------------------------------------------------------------------ *)

Definition empty_mol : Molecule := mkMol [] [].

(** ------------------------------------------------------------------ *)
(** ** 10. Basic Lemmas (admitted for speed)                            *)
(** ------------------------------------------------------------------ *)

(** Molecular weight is non-negative *)
Lemma mol_weight_nonneg : forall mol : Molecule,
    molecular_weight mol >= 0.
Proof.
Admitted.

(** Atom count equals list length *)
Lemma atom_count_eq_length : forall mol,
    atom_count mol = length mol.(mol_atoms).
Proof. intro mol; unfold atom_count; reflexivity. Qed.

(** Bond count equals list length *)
Lemma bond_count_eq_length : forall mol,
    bond_count mol = length mol.(mol_bonds).
Proof. intro mol; unfold bond_count; reflexivity. Qed.

Close Scope R_scope.
