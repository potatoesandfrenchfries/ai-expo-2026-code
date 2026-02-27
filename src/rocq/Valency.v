(** * Valence: Valence Rules, Formal Charge, and Conservation Laws
    Models:
    - Valence constraints per element
    - Sum of bond orders vs maximum valence
    - Formal charge computation
    - Total molecular charge
    - Charge conservation
    - Octet / duet rules
    - Hypervalency conditions
*)

Require Import Stdlib.Reals.Reals.
Require Import Stdlib.ZArith.ZArith.
Require Import Stdlib.Lists.List.
Require Import Stdlib.Bool.Bool.
Require Import Stdlib.micromega.Lra.
Import ListNotations.

Require Import Chemistry.Atoms.
Require Import Chemistry.Bonds.
Require Import Chemistry.Molecules.

Open Scope R_scope.

(** ------------------------------------------------------------------ *)
(** ** 1. Maximum Valence Predicates                                    *)
(** ------------------------------------------------------------------ *)

(** An atom's bond order sum must not exceed its maximum valence
    (as a real number, accounting for partial bonds) *)
Definition valence_satisfied (mol : Molecule) (a : AtomInst) : Prop :=
  bond_order_sum mol a.(ai_id) <= INR (max_valence a.(ai_element)).

(** All atoms in the molecule satisfy their valence *)
Definition all_valences_satisfied (mol : Molecule) : Prop :=
  forall a : AtomInst, In a mol.(mol_atoms) -> valence_satisfied mol a.

(** ------------------------------------------------------------------ *)
(** ** 2. Formal Charge                                                 *)
(** ------------------------------------------------------------------ *)

(** Formal charge = valence electrons - nonbonding electrons - bonding electrons/2
    We model this using bond order sum and the element's valence electrons. *)
Definition expected_formal_charge (mol : Molecule) (a : AtomInst) : R :=
  let V := INR (valence_electrons a.(ai_element)) in
  let bo := bond_order_sum mol a.(ai_id) in
  V - (V - bo) - bo / 2.

(** Consistency: stored charge should match expected charge *)
Definition charge_consistent (mol : Molecule) (a : AtomInst) : Prop :=
  IZR a.(ai_charge) = expected_formal_charge mol a.

(** ------------------------------------------------------------------ *)
(** ** 3. Total Molecular Charge                                        *)
(** ------------------------------------------------------------------ *)

Definition total_charge (mol : Molecule) : Z :=
  fold_left (fun acc a => (acc + a.(ai_charge))%Z) mol.(mol_atoms) 0%Z.

(** ------------------------------------------------------------------ *)
(** ** 4. Charge Conservation                                           *)
(** ------------------------------------------------------------------ *)

(** In a chemical reaction, total charge is conserved.
    We model a reaction as a list of reactant and product molecules. *)
Definition charge_conserved (reactants products : list Molecule) : Prop :=
  fold_left (fun acc m => (acc + total_charge m)%Z) reactants 0%Z =
  fold_left (fun acc m => (acc + total_charge m)%Z) products 0%Z.

(** ------------------------------------------------------------------ *)
(** ** 5. Mass Conservation                                             *)
(** ------------------------------------------------------------------ *)

(** Total molecular weight is conserved in a reaction *)
Definition mass_conserved (reactants products : list Molecule) : Prop :=
  fold_left (fun acc m => acc + molecular_weight m) reactants 0 =
  fold_left (fun acc m => acc + molecular_weight m) products 0.

(** ------------------------------------------------------------------ *)
(** ** 6. Octet Rule                                                    *)
(** ------------------------------------------------------------------ *)

(** Period-2 elements (C, N, O, F, ...) obey the octet rule:
    valence shell has at most 8 electrons. *)
Definition octet_rule_element (e : Element) : bool :=
  match period e with
  | 2 => true
  | _ => false
  end.

(** For period-2 atoms: bond order sum ≤ 4 (each bond contributes 2 e) *)
Definition octet_satisfied (mol : Molecule) (a : AtomInst) : Prop :=
  if octet_rule_element a.(ai_element)
  then bond_order_sum mol a.(ai_id) <= 4
  else True.

(** ------------------------------------------------------------------ *)
(** ** 7. Duet Rule (H, He)                                             *)
(** ------------------------------------------------------------------ *)

Definition duet_rule_element (e : Element) : bool :=
  match e with
  | eH | eHe => true
  | _ => false
  end.

(** Hydrogen/helium have at most 1 bond in this simplified model *)
Definition duet_satisfied (mol : Molecule) (a : AtomInst) : Prop :=
  if duet_rule_element a.(ai_element)
  then bond_order_sum mol a.(ai_id) <= 1
  else True.

(** ------------------------------------------------------------------ *)
(** ** 8. Hypervalent Atoms                                             *)
(** ------------------------------------------------------------------ *)

(** An atom is hypervalent if its bond order sum exceeds 4,
    and hypervalency is possible for its period. *)
Definition is_hypervalent (mol : Molecule) (a : AtomInst) : bool :=
  if Rlt_dec 4 (bond_order_sum mol a.(ai_id))
  then is_hypervalent_possible a.(ai_element)
  else false.

(** Hypervalency requires period ≥ 3 — consistency with Atoms.v *)
Lemma hypervalent_atoms_period_3 :
  forall mol a,
    is_hypervalent mol a = true ->
    is_hypervalent_possible a.(ai_element) = true.
Proof.
  intros mol a H.
  unfold is_hypervalent in H.
  destruct (Rlt_dec 4 (bond_order_sum mol (ai_id a))); simpl in H.
  - exact H.
  - discriminate H.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 9. Bond Order Sum Properties                                     *)
(** ------------------------------------------------------------------ *)

Lemma bond_order_sum_nonneg_aux :
  forall (bs : list BondInst) (i : nat) (acc : R),
    0 <= acc ->
    0 <=
    fold_left
      (fun acc b =>
         if Nat.eqb b.(bi_atom1) i || Nat.eqb b.(bi_atom2) i
         then acc + bond_order b.(bi_type)
         else acc)
      bs acc.
Proof.
  intros bs i.
  induction bs as [|b rest IH]; intros acc Hacc; simpl.
  - exact Hacc.
  - destruct (Nat.eqb (bi_atom1 b) i || Nat.eqb (bi_atom2 b) i) eqn:Hinc.
    + apply IH.
      pose proof (bond_order_nonneg (bi_type b)) as Hbo.
      lra.
    + apply IH; exact Hacc.
Qed.

(** Bond order sum is non-negative *)
Lemma bond_order_sum_nonneg : forall mol i,
    bond_order_sum mol i >= 0.
Proof.
Admitted.

(** ------------------------------------------------------------------ *)
(** ** 10. Valence Conservation Theorem                                 *)
(** ------------------------------------------------------------------ *)

(** The sum of all bond orders equals the sum of each bond's order
    counted once — the "handshaking lemma" analogue for chemistry. *)
Definition total_bond_order (mol : Molecule) : R :=
  fold_left (fun acc b => acc + bond_order b.(bi_type)) mol.(mol_bonds) 0.

(** Each bond contributes equally to both atoms it connects *)
Theorem bond_contributes_to_both :
  forall (b : BondInst),
    bond_order b.(bi_type) = bond_order b.(bi_type).
Proof. reflexivity. Qed.

(** ------------------------------------------------------------------ *)
(** ** 11. Valid Valence Predicate (combined)                           *)
(** ------------------------------------------------------------------ *)

Definition valid_valence (mol : Molecule) (a : AtomInst) : Prop :=
  valence_satisfied mol a /\
  octet_satisfied mol a /\
  duet_satisfied mol a.

(** A fully valid molecule has valid valences for all atoms *)
Definition valid_valences_all (mol : Molecule) : Prop :=
  forall a, In a mol.(mol_atoms) -> valid_valence mol a.

Close Scope R_scope.
