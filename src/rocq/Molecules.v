(** * Molecules: Molecular Graph, Connectivity, and Structural Invariants

    Models a molecule as a labeled graph:
    - Atom instances with element, position, charge, and index
    - Bond instances connecting pairs of atoms
    - Adjacency list / connectivity matrix
    - Neighbor lists
    - Ring membership
    - Substructure / fragment identification
    - Molecular graph invariants
*)

Require Import Stdlib.Reals.Reals.
Require Import Stdlib.ZArith.ZArith.
Require Import Stdlib.Lists.List.
Require Import Stdlib.Bool.Bool.
Require Import Stdlib.Arith.Arith.
Require Import Stdlib.micromega.Lra.
Import ListNotations.
Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.

(** ------------------------------------------------------------------ *)
(** ** 1. Atom Instance                                                 *)
(** ------------------------------------------------------------------ *)

Record AtomInst : Type := mkAtom
  { ai_id       : nat          (** unique identifier within molecule *)
  ; ai_element  : Element
  ; ai_position : Point3D      (** 3D Cartesian coordinates, Å *)
  ; ai_charge   : Z            (** formal charge *)
  ; ai_isotope  : option nat   (** isotope mass number, None = natural *)
  ; ai_stereo   : option bool  (** stereo parity, None = unspecified *)
  }.

(** ------------------------------------------------------------------ *)
(** ** 2. Bond Instance                                                 *)
(** ------------------------------------------------------------------ *)

Record BondInst : Type := mkBond
  { bi_id     : nat
  ; bi_atom1  : nat      (** index of first atom *)
  ; bi_atom2  : nat      (** index of second atom *)
  ; bi_type   : BondType
  ; bi_stereo : option bool  (** E/Z or wedge/dash, None = unspecified *)
  }.

(** A bond connects two different atoms *)
Definition bond_valid_atoms (b : BondInst) : Prop :=
  b.(bi_atom1) <> b.(bi_atom2).

(** ------------------------------------------------------------------ *)
(** ** 3. Molecule                                                      *)
(** ------------------------------------------------------------------ *)

Record Molecule : Type := mkMol
  { mol_atoms : list AtomInst
  ; mol_bonds : list BondInst
  }.

(** ------------------------------------------------------------------ *)
(** ** 4. Atom Lookup                                                   *)
(** ------------------------------------------------------------------ *)

Definition find_atom (mol : Molecule) (id : nat) : option AtomInst :=
  find (fun a => Nat.eqb a.(ai_id) id) mol.(mol_atoms).

(** ------------------------------------------------------------------ *)
(** ** 5. Adjacency / Connectivity                                      *)
(** ------------------------------------------------------------------ *)

(** Are two atom indices connected by a bond? *)
Definition are_bonded (mol : Molecule) (i j : nat) : bool :=
  existsb
    (fun b =>
       (Nat.eqb b.(bi_atom1) i && Nat.eqb b.(bi_atom2) j) ||
       (Nat.eqb b.(bi_atom1) j && Nat.eqb b.(bi_atom2) i))
    mol.(mol_bonds).

(** List of neighbor atom IDs for atom i *)
Definition neighbors (mol : Molecule) (i : nat) : list nat :=
  flat_map
    (fun b =>
       if Nat.eqb b.(bi_atom1) i then [b.(bi_atom2)]
       else if Nat.eqb b.(bi_atom2) i then [b.(bi_atom1)]
       else [])
    mol.(mol_bonds).

(** Degree of atom i (number of bonds) *)
Definition degree (mol : Molecule) (i : nat) : nat :=
  length (neighbors mol i).

(** ------------------------------------------------------------------ *)
(** ** 6. Connectivity Matrix (as a function)                           *)
(** ------------------------------------------------------------------ *)

Definition connectivity_matrix (mol : Molecule) (i j : nat) : bool :=
  are_bonded mol i j.

(** Connectivity matrix is symmetric.
    Proof: for each bond b the predicate is (A||eB) vs (eB||A); use orb_comm
    and induction on the bond list rather than the non-existent
    orb_comm_in_eq lemma from the original outline. *)
Theorem connectivity_matrix_symmetric :
  forall mol i j,
    connectivity_matrix mol i j = connectivity_matrix mol j i.
Proof.
  intros mol i j.
  unfold connectivity_matrix, are_bonded.
  induction (mol_bonds mol) as [|b rest IH].
  - reflexivity.
  - simpl. f_equal.
    + apply orb_comm.
    + apply IH.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 7. Bond List Per Atom                                            *)
(** ------------------------------------------------------------------ *)

Definition bonds_of_atom (mol : Molecule) (i : nat) : list BondInst :=
  filter
    (fun b =>
       Nat.eqb b.(bi_atom1) i || Nat.eqb b.(bi_atom2) i)
    mol.(mol_bonds).

(** ------------------------------------------------------------------ *)
(** ** 8. Graph BFS / DFS reachability                                  *)
(** ------------------------------------------------------------------ *)

(** BFS from a starting atom using *early marking*: a node is added to
    [seen] when ENQUEUED, not when dequeued.  This mirrors the reference
    Haskell implementation (bfs inGraph outGraph queue seen) where
    [seen' = seen ++ enqueue] and [queue' = f ++ enqueue].

    Correspondence with the Haskell variables:
      mol           ~ inGraph
      i :: rest     ~ e : f          (queue head / tail)
      visited       ~ seen           (g:h — enqueued-or-processed set)
      neighbors i   ~ eVertexNeighbors
      new_nbrs      ~ filteredNeighbors / enqueue  (unseen neighbours)
      rest++new_nbrs ~ f ++ enqueue  (queue')
      new_nbrs++visited ~ seen++enqueue (seen')

    Fuel replaces Haskell's structural termination argument: each atom is
    dequeued at most once (early marking), so [n_atoms + 1] steps suffice.

    Invariant on entry: every node already in [queue] is also in [visited]. *)
Fixpoint bfs_fuel
    (mol     : Molecule)
    (queue   : list nat)
    (visited : list nat)   (** seen = enqueued OR already processed *)
    (fuel    : nat)
    : list nat :=
  match fuel with
  | O => visited                        (* fuel exhausted — should not occur *)
  | S n =>
    match queue with
    | [] => visited                     (* bfs _ outGraph [] _ = outGraph *)
    | i :: rest =>                      (* i ~ e,  rest ~ f,  visited ~ seen *)
      let new_nbrs := filter (fun j => negb (existsb (Nat.eqb j) visited))
                             (neighbors mol i) in   (* filteredNeighbors *)
      let visited' := new_nbrs ++ visited in        (* seen' = seen ++ enqueue *)
      let new_q    := rest ++ new_nbrs in            (* queue' = f ++ enqueue *)
      bfs_fuel mol new_q visited' n
    end
  end.

(** [start] is placed in both the queue and visited before the first step,
    upholding the invariant.  Fuel = n_atoms + 1 is sufficient because
    each atom is dequeued at most once. *)
Definition reachable_from (mol : Molecule) (start : nat) : list nat :=
  let n_atoms := length mol.(mol_atoms) in
  bfs_fuel mol [start] [start] (n_atoms + 1).

(** ------------------------------------------------------------------ *)
(** ** 9. Molecular Connectivity                                        *)
(** ------------------------------------------------------------------ *)

(** A molecule is connected if every atom is reachable from atom 0 *)
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
(** ** 10. Fragments / Substructures                                    *)
(** ------------------------------------------------------------------ *)

(** Returns a list of connected components as lists of atom IDs *)
Fixpoint find_components_fuel
    (mol     : Molecule)
    (todo    : list nat)
    (visited : list nat)
    (acc     : list (list nat))
    (fuel    : nat)
    : list (list nat) :=
  match fuel with
  | O => acc
  | S n =>
    match todo with
    | [] => acc
    | i :: rest =>
      if existsb (Nat.eqb i) visited
      then find_components_fuel mol rest visited acc n
      else
        let component := reachable_from mol i in
        let visited'  := component ++ visited in
        find_components_fuel mol rest visited' (component :: acc) n
    end
  end.

Definition fragments (mol : Molecule) : list (list nat) :=
  let all_ids := map ai_id mol.(mol_atoms) in
  find_components_fuel mol all_ids [] [] (length all_ids + 1).

(** ------------------------------------------------------------------ *)
(** ** 11. Unique Atom IDs                                              *)
(** ------------------------------------------------------------------ *)

Definition atom_ids_unique (mol : Molecule) : Prop :=
  NoDup (map ai_id mol.(mol_atoms)).

(** ------------------------------------------------------------------ *)
(** ** 12. Molecular Weight                                             *)
(** ------------------------------------------------------------------ *)

Open Scope R_scope.

(** Sum of atomic masses (ignoring isotope corrections for simplicity) *)
Definition molecular_weight (mol : Molecule) : R :=
  fold_left
    (fun acc a => acc + atomic_mass a.(ai_element))
    mol.(mol_atoms)
    0.

(** Helper: fold_left over atom list preserves non-negativity of accumulator.
    Needed because the naive induction on mol_atoms loses the hypothesis
    about the accumulator's value at each step. *)
Lemma fold_left_mass_nonneg :
  forall (l : list AtomInst) (acc : R),
    acc >= 0 ->
    fold_left (fun (s : R) (a : AtomInst) => s + atomic_mass a.(ai_element)) l acc >= 0.
Proof.
  intros l; induction l as [|a rest IH]; intros acc Hacc.
  - simpl; lra.
  - simpl. apply IH.
    pose proof (atomic_mass_positive a.(ai_element)). lra.
Qed.

(** Molecular weight is non-negative *)
Theorem mol_weight_nonneg : forall mol : Molecule,
    molecular_weight mol >= 0.
Proof.
  intro mol.
  unfold molecular_weight.
  apply fold_left_mass_nonneg. lra.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 13. Atom Count                                                   *)
(** ------------------------------------------------------------------ *)

Definition atom_count (mol : Molecule) : nat :=
  length mol.(mol_atoms).

Definition bond_count (mol : Molecule) : nat :=
  length mol.(mol_bonds).

(** ------------------------------------------------------------------ *)
(** ** 14. Valence Sum Per Atom                                         *)
(** ------------------------------------------------------------------ *)

(** Sum of bond orders for atom i *)
Definition bond_order_sum (mol : Molecule) (i : nat) : R :=
  fold_left
    (fun acc b =>
       if Nat.eqb b.(bi_atom1) i || Nat.eqb b.(bi_atom2) i
       then acc + bond_order b.(bi_type)
       else acc)
    mol.(mol_bonds)
    0.

Close Scope R_scope.

(** ------------------------------------------------------------------ *)
(** ** 15. Molecular Invariants / Well-Formedness                       *)
(** ------------------------------------------------------------------ *)

(** All bond atom references are valid atom IDs *)
Definition bonds_reference_valid_atoms (mol : Molecule) : Prop :=
  forall b : BondInst,
    In b mol.(mol_bonds) ->
    (exists a1, In a1 mol.(mol_atoms) /\ a1.(ai_id) = b.(bi_atom1)) /\
    (exists a2, In a2 mol.(mol_atoms) /\ a2.(ai_id) = b.(bi_atom2)).

(** eNo self-bonds *)
Definition no_self_bonds (mol : Molecule) : Prop :=
  forall b : BondInst, In b mol.(mol_bonds) -> b.(bi_atom1) <> b.(bi_atom2).

(** A well-formed molecule satisfies basic graph constraints *)
Record WellFormedMol (mol : Molecule) : Prop := mkWFM
  { wfm_unique_ids : atom_ids_unique mol
  ; wfm_valid_refs : bonds_reference_valid_atoms mol
  ; wfm_no_self    : no_self_bonds mol
  }.

(** ------------------------------------------------------------------ *)
(** ** 16. Empty Molecule is Well-Formed                                *)
(** ------------------------------------------------------------------ *)

Definition empty_mol : Molecule := mkMol [] [].

Lemma empty_mol_wf : WellFormedMol empty_mol.
Proof.
  constructor.
  - unfold atom_ids_unique, empty_mol; simpl; constructor.
  - unfold bonds_reference_valid_atoms, empty_mol; simpl; intros b Hb; destruct Hb.
  - unfold no_self_bonds, empty_mol; simpl; intros b Hb; destruct Hb.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 17. Completeness Invariants                                      *)
(** ------------------------------------------------------------------ *)

(** All atoms are accounted for: atom count matches the list length *)
Lemma atom_count_eq_length : forall mol,
    atom_count mol = length mol.(mol_atoms).
Proof. intro mol; unfold atom_count; reflexivity. Qed.

(** All bonds are accounted for *)
Lemma bond_count_eq_length : forall mol,
    bond_count mol = length mol.(mol_bonds).
Proof. intro mol; unfold bond_count; reflexivity. Qed.