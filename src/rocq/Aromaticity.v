Require Import Stdlib.Reals.Reals.
Require Import Stdlib.Reals.Rtrigo1.
Require Import Stdlib.Lists.List.
Require Import Stdlib.Arith.Arith.
Require Import Stdlib.micromega.Lra.
Require Import Stdlib.micromega.Lia.
Import ListNotations.

Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.
Require Import Chemistry.Molecules.



(** ------------------------------------------------------------------ *)
(** ** 1. π Electron Count                                              *)
(** ------------------------------------------------------------------ *)

(** Number of π electrons contributed by an element in an aromatic ring
    (simplified model):
    - sp² C, N (pyridine-like): contributes 1 π electron
    - N with lone pair (pyrrole-like), O, S in ring: contributes 2 π electrons
    - default: 0 *)
Inductive PiContribution : Type :=
  | PiZero
  | PiOne
  | PiTwo.

Definition pi_contribution_val (pc : PiContribution) : nat :=
  match pc with
  | PiZero => 0
  | PiOne  => 1
  | PiTwo  => 2
  end.

(** ------------------------------------------------------------------ *)
(** ** 2. Hückel's Rule                                                 *)
(** ------------------------------------------------------------------ *)

(** A count of π electrons satisfies Hückel's rule iff it equals 4n+2
    for some non-negative integer n. *)
Definition huckel_4n2 (pi_count : nat) : Prop :=
  exists n : nat, pi_count = 4 * n + 2.

(** Hückel's rule is decidable *)
Definition huckel_4n2_dec (pi_count : nat) : bool :=
  match pi_count mod 4 with
  | 2 => true
  | _ => false
  end.

Theorem huckel_dec_correct : forall n : nat,
    huckel_4n2 n <-> huckel_4n2_dec n = true.
Proof.
  intro n; split.
  - intros [k Hk]. subst.
    unfold huckel_4n2_dec.
    rewrite Nat.add_mod by lia.
    rewrite Nat.mul_mod by lia.
    rewrite Nat.mod_same by lia.
    rewrite Nat.mod_0_l by lia.
    simpl; reflexivity.
  - intro H.
    unfold huckel_4n2_dec in H.
    destruct (n mod 4) as [| [| [|r]]] eqn:Hmod; simpl in H; try discriminate.
    exists (n / 4).
    rewrite (Nat.div_mod n 4) at 1 by lia.
    rewrite Hmod.
    lia.
Qed.

(** Benzene has 6 π electrons: satisfies Hückel's rule (n=1) *)
Lemma benzene_huckel : huckel_4n2 6.
Proof. exists 1; reflexivity. Qed.

(** Naphthalene has 10 π electrons: satisfies Hückel's rule (n=2) *)
Lemma naphthalene_huckel : huckel_4n2 10.
Proof. exists 2; reflexivity. Qed.

(** Anthracene has 14 π electrons: satisfies Hückel's rule (n=3) *)
Lemma anthracene_huckel : huckel_4n2 14.
Proof. exists 3; reflexivity. Qed.

(** Cyclopentadienyl anion: 6 π electrons, Hückel satisfied *)
Lemma cyclopentadienyl_anion_huckel : huckel_4n2 6.
Proof. exists 1; reflexivity. Qed.

(** Tropylium cation: 6 π electrons, Hückel satisfied *)
Lemma tropylium_huckel : huckel_4n2 6.
Proof. exists 1; reflexivity. Qed.


(** ------------------------------------------------------------------ *)
(** ** 3. Anti-aromaticity (4n π electrons)                             *)
(** ------------------------------------------------------------------ *)

Definition anti_aromatic_4n (pi_count : nat) : Prop :=
  exists n : nat, n > 0 /\ pi_count = 4 * n.

Definition anti_aromatic_4n_dec (pi_count : nat) : bool :=
  match pi_count with
  | 0 => false
  | _ =>
    match pi_count mod 4 with
    | 0 => true
    | _ => false
    end
  end.

(** Cyclobutadiene: 4 π electrons → anti-aromatic *)
Lemma cyclobutadiene_antiaromatic : anti_aromatic_4n 4.
Proof. exists 1; split; [lia | reflexivity]. Qed.

(** Hückel and anti-aromatic are mutually exclusive *)
Theorem huckel_not_antiaromatic : forall n : nat,
    huckel_4n2 n -> ~ anti_aromatic_4n n.
Proof.
  intros n [k Hk] [m [Hm Hm']].
  subst.
  lia.
Qed.


(** ------------------------------------------------------------------ *)
(** ** 4. Planarity Requirement                                         *)
(** ------------------------------------------------------------------ *)
Open Scope R_scope.

(** A ring is planar if all its atoms lie within ε of a common plane.
    We define planarity as maximum z-deviation after centering ≤ threshold. *)
Definition max_z_deviation (positions : list Point3D) : R :=
  match positions with
  | [] => 0
  | first :: rest =>
    let z_avg :=
      fold_left (fun acc p => acc + p.z) positions 0 /
      INR (length positions)
    in
    fold_left
      (fun acc p => Rmax acc (Rabs (p.z - z_avg)))
      positions
      0
  end.

(** Aromaticity requires near-planarity: deviation ≤ 0.1 Å *)
Definition planarity_threshold : R := 0.1.

Definition is_planar (positions : list Point3D) : bool :=
  if Rle_dec (max_z_deviation positions) planarity_threshold then true else false.


(** ------------------------------------------------------------------ *)
(** ** 5. Bond Length Equalization                                      *)
(** ------------------------------------------------------------------ *)

(** In an aromatic ring, C-C bonds should be equalized (Å) *)
Definition aromatic_CC_min : R := 1.38.
Definition aromatic_CC_max : R := 1.42.

Definition is_aromatic_bond_length (len : R) : bool :=
  if Rle_dec aromatic_CC_min len then
    if Rle_dec len aromatic_CC_max then true else false
  else false.
(** ------------------------------------------------------------------ *)
(** ** 6. Full Aromaticity Predicate                                    *)
(** ------------------------------------------------------------------ *)

(** A ring is aromatic if:
    1. π electron count satisfies Hückel's rule
    2. The ring is planar (positions lie in a plane)
    3. Bond lengths are equalized (for C-C, checked separately) *)
Record AromaticRing : Type := mkAromaticRing
  { ar_atom_ids  : list nat       (** ordered ring atom IDs *)
  ; ar_pi_count  : nat            (** total π electrons *)
  ; ar_positions : list Point3D   (** 3D positions of ring atoms *)
  }.

Definition is_aromatic_ring (ar : AromaticRing) : bool :=
  huckel_4n2_dec ar.(ar_pi_count) &&
  is_planar ar.(ar_positions).

(** ------------------------------------------------------------------ *)
(** ** 7. Specific Aromatic Systems                                     *)
(** ------------------------------------------------------------------ *)

(** Benzene: 6-membered ring with 6 π electrons *)
Definition benzene_pi_count  : nat := 6.
Definition benzene_ring_size : nat := 6.

(** Pyridine: 6-membered ring, N contributes 1 π electron *)
Definition pyridine_pi_count : nat := 6.

(** Pyrrole: 5-membered ring, NH contributes 2 π electrons *)
Definition pyrrole_pi_count  : nat := 6.

(** Furan: 5-membered ring, O contributes 2 π electrons *)
Definition furan_pi_count    : nat := 6.

(** Thiophene: 5-membered ring, S contributes 2 π electrons *)
Definition thiophene_pi_count : nat := 6.

Close Scope R_scope.
(** All these satisfy Hückel's rule *)
Lemma benzene_aromatic   : huckel_4n2 benzene_pi_count.   Proof. exists 1; reflexivity. Qed.
Lemma pyridine_aromatic  : huckel_4n2 pyridine_pi_count.  Proof. exists 1; reflexivity. Qed.
Lemma pyrrole_aromatic   : huckel_4n2 pyrrole_pi_count.   Proof. exists 1; reflexivity. Qed.
Lemma furan_aromatic     : huckel_4n2 furan_pi_count.     Proof. exists 1; reflexivity. Qed.
Lemma thiophene_aromatic : huckel_4n2 thiophene_pi_count. Proof. exists 1; reflexivity. Qed.

(** ------------------------------------------------------------------ *)
(** ** 8. Aromatic Stabilization Energy                                 *)
(** ------------------------------------------------------------------ *)

(** Benzene resonance stabilization energy ≈ 150 kJ/mol *)
Definition benzene_resonance_energy : R := 150.0.

(** ------------------------------------------------------------------ *)
(** ** 9. Magnetic Ring Current (NMR chemical shifts)                   *)
(** ------------------------------------------------------------------ *)

(** Aromatic protons experience deshielding (δ ≈ 7-9 ppm) due to ring current *)
Definition aromatic_proton_shift_min : R := 7.0.
Definition aromatic_proton_shift_max : R := 9.0.

(** Antiaromatic protons experience shielding (can be negative δ) *)
Definition antiaromatic_proton_shift_max : R := 5.0.

