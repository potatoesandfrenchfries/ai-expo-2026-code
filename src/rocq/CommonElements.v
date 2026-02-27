(** * CommonElements: Optimized subset of elements for drug-like molecules
    
    This module provides a minimal set of elements commonly found in organic
    and medicinal chemistry, dramatically reducing compilation time compared
    to the full 118-element periodic table.
    
    Includes: H, C, N, O, S, P, F, Cl, Br, I (10 elements)
*)

Require Import Stdlib.Reals.Reals.
Require Import Stdlib.ZArith.ZArith.
Require Import Stdlib.micromega.Lia.
Require Import Stdlib.micromega.Lra.
Require Import Stdlib.Lists.List.
Import ListNotations.

Open Scope R_scope.

(** ** 1. Common Element Types (10 elements for drug-like molecules) *)

Inductive CommonElement : Type :=
  | eH   (* Hydrogen *)
  | eC   (* Carbon *)
  | eN   (* Nitrogen *)
  | eO   (* Oxygen *)
  | eS   (* Sulfur *)
  | eP   (* Phosphorus *)
  | eF   (* Fluorine *)
  | eCl  (* Chlorine *)
  | eBr  (* Bromine *)
  | eI.  (* Iodine *)

(** Element equality is decidable *)
Lemma common_element_eq_dec : forall (a b : CommonElement), {a = b} + {a <> b}.
Proof. decide equality. Defined.

(** ** 2. Atomic Number *)

Definition atomic_number (e : CommonElement) : nat :=
  match e with
  | eH  => 1
  | eC  => 6
  | eN  => 7
  | eO  => 8
  | eF  => 9
  | eP  => 15
  | eS  => 16
  | eCl => 17
  | eBr => 35
  | eI  => 53
  end.

(** ** 3. Atomic Mass (in Daltons, u) *)

Definition atomic_mass (e : CommonElement) : R :=
  match e with
  | eH  => 1.008
  | eC  => 12.011
  | eN  => 14.007
  | eO  => 15.999
  | eF  => 18.998
  | eP  => 30.974
  | eS  => 32.065
  | eCl => 35.453
  | eBr => 79.904
  | eI  => 126.90
  end.

(** ** 4. Van der Waals Radius (in Angstroms, Å) *)

Definition van_der_waals_radius (e : CommonElement) : R :=
  match e with
  | eH  => 1.20
  | eC  => 1.70
  | eN  => 1.55
  | eO  => 1.52
  | eF  => 1.47
  | eP  => 1.80
  | eS  => 1.80
  | eCl => 1.75
  | eBr => 1.85
  | eI  => 1.98
  end.

(** ** 5. Covalent Radius (in Angstroms, Å) *)

Definition covalent_radius (e : CommonElement) : R :=
  match e with
  | eH  => 0.31
  | eC  => 0.77
  | eN  => 0.71
  | eO  => 0.66
  | eF  => 0.57
  | eP  => 1.07
  | eS  => 1.05
  | eCl => 1.02
  | eBr => 1.20
  | eI  => 1.39
  end.

(** ** 6. Pauling Electronegativity *)

Definition electronegativity (e : CommonElement) : R :=
  match e with
  | eH  => 2.20
  | eC  => 2.55
  | eN  => 3.04
  | eO  => 3.44
  | eF  => 3.98
  | eP  => 2.19
  | eS  => 2.58
  | eCl => 3.16
  | eBr => 2.96
  | eI  => 2.66
  end.

(** ** 7. First Ionization Energy (kJ/mol) *)

Definition ionization_energy (e : CommonElement) : R :=
  match e with
  | eH  => 1312.0
  | eC  => 1086.5
  | eN  => 1402.3
  | eO  => 1313.9
  | eF  => 1681.0
  | eP  => 1011.8
  | eS  => 999.6
  | eCl => 1251.2
  | eBr => 1139.9
  | eI  => 1008.4
  end.

(** ** 8. Electron Affinity (kJ/mol, positive = exothermic) *)

Definition electron_affinity (e : CommonElement) : R :=
  match e with
  | eH  => 72.8
  | eC  => 121.8
  | eN  => 0.0
  | eO  => 141.0
  | eF  => 328.2
  | eP  => 72.0
  | eS  => 200.4
  | eCl => 348.6
  | eBr => 324.6
  | eI  => 295.2
  end.

(** ** 9. Period (row in periodic table) *)

Definition period (e : CommonElement) : nat :=
  match e with
  | eH  => 1
  | eC | eN | eO | eF => 2
  | eP | eS | eCl => 3
  | eBr => 4
  | eI  => 5
  end.

(** ** 10. Group (column in periodic table, IUPAC 1-18) *)

Definition group (e : CommonElement) : nat :=
  match e with
  | eH  => 1
  | eC  => 14
  | eN  => 15
  | eO  => 16
  | eF | eCl | eBr | eI => 17
  | eP  => 15
  | eS  => 16
  end.

(** ** 11. Valence Electrons *)

Definition valence_electrons (e : CommonElement) : nat :=
  match e with
  | eH  => 1
  | eC  => 4
  | eN  => 5
  | eO  => 6
  | eF  => 7
  | eP  => 5
  | eS  => 6
  | eCl => 7
  | eBr => 7
  | eI  => 7
  end.

(** ** 12. Maximum Valence *)

Definition max_valence (e : CommonElement) : nat :=
  match e with
  | eH  => 1
  | eC  => 4
  | eN  => 3
  | eO  => 2
  | eF  => 1
  | eP  => 5
  | eS  => 6
  | eCl => 7
  | eBr => 5
  | eI  => 7
  end.

(** ** 13. Hypervalent atoms (period ≥ 3 allowing expanded octet) *)

Definition is_hypervalent_possible (e : CommonElement) : bool :=
  match period e with
  | 0 | 1 | 2 => false
  | _          => true
  end.

(** ** 14. Basic Properties/Lemmas *)

(** Atomic number is positive for all elements *)
Lemma atomic_number_positive : forall e : CommonElement, (atomic_number e >= 1)%nat.
Proof.
  intro e; destruct e; simpl; lia.
Qed.

(** Atomic mass is positive *)
Lemma atomic_mass_positive : forall e : CommonElement, atomic_mass e > 0.
Proof.
  intro e; destruct e; simpl; lra.
Qed.

(** Van der Waals radius is positive *)
Lemma vdw_radius_positive : forall e : CommonElement, van_der_waals_radius e > 0.
Proof.
  intro e; destruct e; simpl; lra.
Qed.

(** Electronegativity is non-negative *)
Lemma electronegativity_nonneg : forall e : CommonElement, electronegativity e >= 0.
Proof.
  intro e; destruct e; simpl; lra.
Qed.

(** Fluorine has the highest electronegativity (Pauling scale) *)
Lemma fluorine_max_electronegativity :
  forall e : CommonElement, electronegativity e <= electronegativity eF.
Proof.
  intro e; destruct e; simpl; lra.
Qed.

(** Period is between 1 and 5 for common elements *)
Lemma period_bounds : forall e : CommonElement, (1 <= period e /\ period e <= 5)%nat.
Proof.
  intro e; destruct e; simpl; lia.
Qed.

Close Scope R_scope.
