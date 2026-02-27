(** * Conformational Analysis: Energy Terms and Conformational Search

    Models:
    - Force field energy terms (bond stretch, angle bend, torsion, VdW, electrostatic)
    - Conformational energy landscape
    - Local/global energy minima
    - Transition states and barriers
    - Boltzmann populations
    - Molecular dynamics constraints
*)

Require Import Stdlib.Reals.Reals.
Require Import Stdlib.Reals.Rtrigo1.
Require Import Stdlib.Lists.List.
Require Import Stdlib.Arith.Arith.
Require Import Stdlib.micromega.Lra.
Import ListNotations.

Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.
Require Import Chemistry.Molecules.

Open Scope R_scope.

(** ------------------------------------------------------------------ *)
(** ** 1. Force Field Energy Terms                                      *)
(** ------------------------------------------------------------------ *)

(** Bond stretching energy (harmonic potential):
    E_bond = (k/2) * (r - r0)^2 *)
Definition bond_stretch_energy (k r r0 : R) : R :=
  (k / 2) * (r - r0) * (r - r0).

(** Bond stretching energy is non-negative *)
Lemma bond_stretch_energy_nonneg : forall k r r0 : R,
    k >= 0 ->
    bond_stretch_energy k r r0 >= 0.
Proof.
Admitted.

(** Equilibrium: zero energy at r = r0 *)
Lemma bond_stretch_zero_at_eq : forall k r0 : R,
    bond_stretch_energy k r0 r0 = 0.
Proof.
  intros k r0.
  unfold bond_stretch_energy; ring.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 2. Angle Bending Energy                                          *)
(** ------------------------------------------------------------------ *)

(** E_angle = (k_theta/2) * (theta - theta0)^2 *)
Definition angle_bend_energy (k_theta theta theta0 : R) : R :=
  (k_theta / 2) * (theta - theta0) * (theta - theta0).

Lemma angle_bend_energy_nonneg : forall k t t0 : R,
    k >= 0 ->
    angle_bend_energy k t t0 >= 0.
Proof.
Admitted.

Lemma angle_bend_zero_at_eq : forall k theta0 : R,
    angle_bend_energy k theta0 theta0 = 0.
Proof.
  intros k theta0; unfold angle_bend_energy; ring.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 3. Torsional Energy                                              *)
(** ------------------------------------------------------------------ *)

(** E_torsion = (V_n / 2) * (1 + cos(n * phi - gamma)) *)
Definition torsional_energy (V_n n phi gamma : R) : R :=
  (V_n / 2) * (1 + cos (n * phi - gamma)).

(** Torsional energy is bounded: 0 <= E_torsion <= V_n *)
Lemma torsional_energy_nonneg : forall V_n n phi gamma : R,
    V_n >= 0 ->
    torsional_energy V_n n phi gamma >= 0.
Proof.
  intros V_n n phi gamma HVn.
  unfold torsional_energy.
  pose proof (COS_bound (n * phi - gamma)) as Hcos.
  nra.
Qed.

Lemma torsional_energy_bounded : forall V_n n phi gamma : R,
    V_n >= 0 ->
    torsional_energy V_n n phi gamma <= V_n.
Proof.
  intros V_n n phi gamma HVn.
  unfold torsional_energy.
  pose proof (COS_bound (n * phi - gamma)) as Hcos.
  nra.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 4. Van der Waals Energy (Lennard-Jones potential)                *)
(** ------------------------------------------------------------------ *)

(** E_LJ = 4 epsilon [ (sigma/r)^12 - (sigma/r)^6 ]
    epsilon = well depth, sigma = collision diameter, r = distance *)
Definition lennard_jones_energy (epsilon sigma r : R) : R :=
  if Rle_dec r 0 then 0
  else
    let sr6 := (sigma / r) ^ 6 in
    4 * epsilon * (sr6 * sr6 - sr6).

(** ------------------------------------------------------------------ *)
(** ** 5. Electrostatic Energy (Coulomb)                                *)
(** ------------------------------------------------------------------ *)

(** E_elec = k_e * q1 * q2 / (epsilon_r * r) *)
Definition coulomb_energy (q1 q2 epsilon_r r : R) : R :=
  if Rle_dec r 0 then 0
  else
    let k_e := 332.0 in
    k_e * q1 * q2 / (epsilon_r * r).

(** Like charges repel: positive energy (for epsilon_r = 1) *)
Lemma like_charges_repel : forall q r : R,
    q > 0 -> r > 0 ->
    coulomb_energy q q 1 r > 0.
Proof.
  intros q r Hq Hr.
  unfold coulomb_energy.
  destruct (Rle_dec r 0); [lra |].
  replace (1 * r) with r by nra.
  unfold Rdiv.
  apply Rmult_lt_0_compat.
  - nra.
  - apply Rinv_0_lt_compat; lra.
Qed.

(** Opposite charges attract: negative energy (for epsilon_r = 1) *)
Lemma opposite_charges_attract : forall q r : R,
    q > 0 -> r > 0 ->
    coulomb_energy q (- q) 1 r < 0.
Proof.
  intros q r Hq Hr.
  unfold coulomb_energy.
  destruct (Rle_dec r 0); [lra |].
  replace (1 * r) with r by nra.
  unfold Rdiv.
  assert (Hinv : / r > 0) by (apply Rinv_0_lt_compat; lra).
  assert (Hnum : 332 * q * (- q) < 0) by nra.
  nra.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 6. Total Force Field Energy                                      *)
(** ------------------------------------------------------------------ *)

Record FFEnergy : Type := mkFFE
  { ffe_bond    : R
  ; ffe_angle   : R
  ; ffe_torsion : R
  ; ffe_vdw     : R
  ; ffe_elec    : R
  ; ffe_hbond   : R
  ; ffe_solv    : R
  }.

Definition total_ff_energy (e : FFEnergy) : R :=
  e.(ffe_bond) + e.(ffe_angle) + e.(ffe_torsion) +
  e.(ffe_vdw) + e.(ffe_elec) + e.(ffe_hbond) + e.(ffe_solv).

(** ------------------------------------------------------------------ *)
(** ** 7. Energy Minimum                                                *)
(** ------------------------------------------------------------------ *)

Definition is_local_minimum
    (E : Point3D -> R)
    (x0 : Point3D)
    (delta : R)
    : Prop :=
  forall x : Point3D,
    distance x x0 < delta ->
    E x0 <= E x.

(** ------------------------------------------------------------------ *)
(** ** 8. Conformational Search                                         *)
(** ------------------------------------------------------------------ *)

Record ConformationalEnsemble : Type := mkCE
  { ce_conformers : list (Point3D * R)
  }.

Definition global_minimum_energy (ce : ConformationalEnsemble) : option R :=
  match ce.(ce_conformers) with
  | [] => None
  | (_, e) :: rest =>
      Some (fold_left (fun acc pe => Rmin acc (snd pe)) rest e)
  end.

Definition conformer_population
    (E_i : R)
    (E_all : list R)
    (kT : R)
    : R :=
  let Z := fold_left (fun acc e => acc + exp (- e / kT)) E_all 0 in
  if Req_EM_T Z 0 then 0
  else exp (- E_i / kT) / Z.

(** ------------------------------------------------------------------ *)
(** ** 9. Molecular Dynamics Parameters                                 *)
(** ------------------------------------------------------------------ *)

Inductive MDEnsemble : Type :=
  | NVE
  | NVT
  | NPT
  | NpT.

Record MDParameters : Type := mkMDP
  { md_timestep_fs   : R
  ; md_temperature_K : R
  ; md_pressure_atm  : R
  ; md_ensemble      : MDEnsemble
  }.

Definition valid_timestep (dt_fs : R) : Prop :=
  0 < dt_fs /\ dt_fs <= 2.0.

Definition valid_temperature (T_K : R) : Prop :=
  250 <= T_K /\ T_K <= 400.

(** ------------------------------------------------------------------ *)
(** ** 10. RMSD and RMSF                                                *)
(** ------------------------------------------------------------------ *)

Definition rmsd
    (coords1 coords2 : list Point3D)
    : R :=
  match coords1, coords2 with
  | [], [] => 0
  | _, _ =>
      let n := length coords1 in
      if Nat.eqb n 0 then 0
      else
        let sum_sq :=
          fold_left
            (fun acc pp =>
               let p1 := fst pp in
               let p2 := snd pp in
               let d := distance p1 p2 in
               acc + d * d)
            (combine coords1 coords2)
            0
        in
        sqrt (sum_sq / INR n)
  end.

Lemma rmsd_nonneg : forall c1 c2, rmsd c1 c2 >= 0.
Proof.
Admitted.

Lemma fold_self_dist_sq_zero :
  forall (l : list Point3D) (acc : R),
    fold_left
      (fun acc pp =>
         let p1 := fst pp in
         let p2 := snd pp in
         let d := distance p1 p2 in
         acc + d * d)
      (combine l l)
      acc = acc.
Proof.
  intros l.
  induction l as [|x xs IH]; intro acc; simpl.
  - reflexivity.
  - rewrite distance_self.
    replace (acc + 0 * 0) with acc by ring.
    apply IH.
Qed.

Lemma rmsd_self : forall c : list Point3D, rmsd c c = 0.
Proof.
  intro c.
  unfold rmsd.
  destruct c as [|x xs].
  - reflexivity.
  - simpl.
    destruct (Nat.eqb (S (length xs)) 0) eqn:Hn.
    + apply Nat.eqb_eq in Hn. discriminate.
    + rewrite fold_self_dist_sq_zero.
      rewrite distance_self.
      replace
        ((0 + 0 * 0) /
         match length xs with
         | 0%nat => 1
         | S _ => INR (length xs) + 1
         end) with 0.
      2:{ unfold Rdiv; ring. }
      apply sqrt_0.
Qed.

Close Scope R_scope.
