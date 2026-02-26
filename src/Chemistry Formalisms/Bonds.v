(** * Bonds: Bond Types, Orders, and Geometric Constraints

    Models chemical bonds including:
    - Bond types (single, double, triple, aromatic, etc.)
    - Bond orders
    - Bond length constraints (per element pair)
    - Bond angle constraints (per hybridization)
    - Dihedral angle constraints
    - Bond polarity and energy
*)

Require Import Stdlib.Reals.Reals.
Require Import Stdlib.Reals.Rtrigo_def.
Require Import Stdlib.Reals.RIneq.
Require Import Stdlib.Lists.List.
Require Import Stdlib.micromega.Lra.
Import ListNotations.
Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.

Open Scope R_scope.

(** ------------------------------------------------------------------ *)
(** ** 1. Bond Types                                                    *)
(** ------------------------------------------------------------------ *)

Inductive BondType : Type :=
  | SingleBond
  | DoubleBond
  | TripleBond
  | AromaticBond
  | CoordinateBond      (** dative / coordinate covalent *)
  | IonicBond
  | HydrogenBond
  | MetallicBond.

(** Bond type equality is decidable *)
Lemma bondtype_eq_dec : forall (a b : BondType), {a = b} + {a <> b}.
Proof. decide equality. Defined.

(** ------------------------------------------------------------------ *)
(** ** 2. Bond Order                                                    *)
(** ------------------------------------------------------------------ *)

(** Bond order as a real number (1, 1.5, 2, 3, …) *)
Definition bond_order (bt : BondType) : R :=
  match bt with
  | SingleBond     => 1
  | DoubleBond     => 2
  | TripleBond     => 3
  | AromaticBond   => 1.5
  | CoordinateBond => 1
  | IonicBond      => 0   (* approximation; ionic has no shared electrons *)
  | HydrogenBond   => 0
  | MetallicBond   => 0
  end.

(** Bond order is non-negative *)
Lemma bond_order_nonneg : forall bt, bond_order bt >= 0.
Proof. intro bt; destruct bt; simpl; lra. Qed.

(** Covalent bonds have positive bond order *)
Definition is_covalent (bt : BondType) : bool :=
  match bt with
  | SingleBond | DoubleBond | TripleBond | AromaticBond | CoordinateBond => true
  | _ => false
  end.

(** ------------------------------------------------------------------ *)
(** ** 3. Bond Length Constraints (in Angstroms)                        *)
(** ------------------------------------------------------------------ *)

(** Each constraint is a (min, max) interval in Å *)
Record BondLengthRange : Type := mkBLR
  { bl_min : R ; bl_max : R }.

(** The bond is valid if its length falls within [bl_min, bl_max] *)
Definition valid_bond_length (blr : BondLengthRange) (len : R) : Prop :=
  blr.(bl_min) <= len /\ len <= blr.(bl_max).

(** Validity is decidable given a concrete length.
    Note: Rle_dec returns sumbool, so we use nested if-then-else. *)
Definition bond_length_in_range (blr : BondLengthRange) (len : R) : bool :=
  if Rle_dec blr.(bl_min) len then
    if Rle_dec len blr.(bl_max) then true
    else false
  else false.

(** Standard bond length ranges (Å) for common element pairs *)

(** C-C bonds *)
Definition bl_CC_single  : BondLengthRange := mkBLR 1.47 1.61.
Definition bl_CC_double  : BondLengthRange := mkBLR 1.30 1.40.
Definition bl_CC_triple  : BondLengthRange := mkBLR 1.18 1.24.
Definition bl_CC_aromatic: BondLengthRange := mkBLR 1.38 1.42.

(** C-N bonds *)
Definition bl_CN_single  : BondLengthRange := mkBLR 1.43 1.50.
Definition bl_CN_double  : BondLengthRange := mkBLR 1.27 1.35.
Definition bl_CN_triple  : BondLengthRange := mkBLR 1.13 1.18.

(** C-O bonds *)
Definition bl_CO_single  : BondLengthRange := mkBLR 1.39 1.46.
Definition bl_CO_double  : BondLengthRange := mkBLR 1.18 1.24.

(** C-H, N-H, O-H *)
Definition bl_CH         : BondLengthRange := mkBLR 1.06 1.12.
Definition bl_NH         : BondLengthRange := mkBLR 0.99 1.05.
Definition bl_OH         : BondLengthRange := mkBLR 0.94 1.00.

(** N-N, N-O, O-O *)
Definition bl_NN_single  : BondLengthRange := mkBLR 1.40 1.50.
Definition bl_NN_double  : BondLengthRange := mkBLR 1.20 1.30.
Definition bl_NO_single  : BondLengthRange := mkBLR 1.36 1.44.
Definition bl_NO_double  : BondLengthRange := mkBLR 1.15 1.25.
Definition bl_OO_single  : BondLengthRange := mkBLR 1.40 1.50.

(** Halogens to C *)
Definition bl_CF         : BondLengthRange := mkBLR 1.30 1.40.
Definition bl_CCl        : BondLengthRange := mkBLR 1.75 1.85.
Definition bl_CBr        : BondLengthRange := mkBLR 1.90 2.00.
Definition bl_CI         : BondLengthRange := mkBLR 2.10 2.20.

(** Sulfur and phosphorus *)
Definition bl_SS         : BondLengthRange := mkBLR 2.00 2.10.
Definition bl_CS         : BondLengthRange := mkBLR 1.75 1.85.
Definition bl_CP         : BondLengthRange := mkBLR 1.80 1.90.
Definition bl_PO         : BondLengthRange := mkBLR 1.55 1.65.
Definition bl_SO         : BondLengthRange := mkBLR 1.45 1.55.

(** The range is well-formed (min ≤ max) for all standard definitions *)
Lemma bl_CC_single_wf   : bl_CC_single.(bl_min) <= bl_CC_single.(bl_max).   Proof. simpl; lra. Qed.
Lemma bl_CC_double_wf   : bl_CC_double.(bl_min) <= bl_CC_double.(bl_max).   Proof. simpl; lra. Qed.
Lemma bl_CC_triple_wf   : bl_CC_triple.(bl_min) <= bl_CC_triple.(bl_max).   Proof. simpl; lra. Qed.
Lemma bl_CC_aromatic_wf : bl_CC_aromatic.(bl_min) <= bl_CC_aromatic.(bl_max). Proof. simpl; lra. Qed.

(** Bond length decreases with increasing bond order for C-C *)
Theorem CC_bond_length_order :
  bl_CC_triple.(bl_max) < bl_CC_double.(bl_min) /\
  bl_CC_double.(bl_max) < bl_CC_single.(bl_min).
Proof.
  unfold bl_CC_triple, bl_CC_double, bl_CC_single; simpl; lra.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 4. Hybridization and Bond Angle Constraints                      *)
(** ------------------------------------------------------------------ *)

Inductive Hybridization : Type :=
  | Sp3   (** tetrahedral, 109.5° *)
  | Sp2   (** trigonal planar, 120° *)
  | Sp    (** linear, 180° *)
  | Other.

(** Ideal bond angle in radians for each hybridization *)
Definition ideal_angle (h : Hybridization) : R :=
  match h with
  | Sp3   => 109.5 * PI / 180   (* ≈ 1.911 rad *)
  | Sp2   => 120.0 * PI / 180   (* ≈ 2.094 rad *)
  | Sp    => PI                  (* 180° *)
  | Other => 0
  end.

(** Tolerance in radians for each hybridization *)
Definition angle_tolerance (h : Hybridization) : R :=
  match h with
  | Sp3   => 5.0 * PI / 180   (* ±5° *)
  | Sp2   => 5.0 * PI / 180
  | Sp    => 3.0 * PI / 180
  | Other => 0
  end.

Record BondAngleRange : Type := mkBAR
  { ba_min : R ; ba_max : R }.

Definition hybridization_angle_range (h : Hybridization) : BondAngleRange :=
  let mid := ideal_angle h in
  let tol := angle_tolerance h in
  mkBAR (mid - tol) (mid + tol).

Definition valid_bond_angle (h : Hybridization) (theta : R) : Prop :=
  let r := hybridization_angle_range h in
  r.(ba_min) <= theta /\ theta <= r.(ba_max).

(** Specific bond angle constraints (in degrees, for readability) *)

(** sp³ tetrahedral: 109.5° ± 5° *)
Definition sp3_angle_min_deg : R := 104.5.
Definition sp3_angle_max_deg : R := 114.5.

(** sp² trigonal planar: 120° ± 5° *)
Definition sp2_angle_min_deg : R := 115.0.
Definition sp2_angle_max_deg : R := 125.0.

(** sp linear: 180° ± 3° *)
Definition sp_angle_min_deg  : R := 177.0.
Definition sp_angle_max_deg  : R := 183.0.

(** Water H-O-H: 104.5° ± 2° *)
Definition water_angle_min_deg : R := 102.5.
Definition water_angle_max_deg : R := 106.5.

(** Ammonia H-N-H: 107° ± 2° *)
Definition ammonia_angle_min_deg : R := 105.0.
Definition ammonia_angle_max_deg : R := 109.0.

(** Ring angles *)
Definition cyclopropane_angle_min : R := 58.0.
Definition cyclopropane_angle_max : R := 62.0.
Definition cyclobutane_angle_min  : R := 85.0.
Definition cyclobutane_angle_max  : R := 95.0.
Definition cyclopentane_angle_min : R := 105.0.
Definition cyclopentane_angle_max : R := 111.0.
Definition cyclohexane_angle_min  : R := 109.0.
Definition cyclohexane_angle_max  : R := 113.0.
Definition benzene_angle_min      : R := 119.0.
Definition benzene_angle_max      : R := 121.0.

(** ------------------------------------------------------------------ *)
(** ** 5. Dihedral Angle Constraints                                    *)
(** ------------------------------------------------------------------ *)

(** Dihedral angles for common conformations (in degrees) *)

(** Staggered ethane: 60°, 180°, 300° *)
Definition ethane_staggered_1 : R :=  60.0.
Definition ethane_staggered_2 : R := 180.0.
Definition ethane_staggered_3 : R := 300.0.

(** Eclipsed ethane: 0°, 120°, 240° *)
Definition ethane_eclipsed_1  : R :=   0.0.
Definition ethane_eclipsed_2  : R := 120.0.
Definition ethane_eclipsed_3  : R := 240.0.

(** Double bond geometry (in degrees) with 5° tolerance *)
Definition trans_double_bond_min : R := 175.0.
Definition trans_double_bond_max : R := 185.0.
Definition cis_double_bond_min   : R := -5.0.
Definition cis_double_bond_max   : R :=  5.0.

(** Peptide bond (trans predominates in proteins) *)
Definition peptide_trans_min : R := 170.0.
Definition peptide_trans_max : R := 190.0.  (* wraps; treat mod 360 *)
Definition peptide_cis_min   : R := -10.0.
Definition peptide_cis_max   : R :=  10.0.

(** ------------------------------------------------------------------ *)
(** ** 6. Bond Energy (approximate, kJ/mol)                             *)
(** ------------------------------------------------------------------ *)

Definition bond_energy_CC (bt : BondType) : R :=
  match bt with
  | SingleBond   => 347.0
  | DoubleBond   => 614.0
  | TripleBond   => 839.0
  | AromaticBond => 507.0   (* ~average of single + double *)
  | _            => 0.0
  end.

Definition bond_energy_CN (bt : BondType) : R :=
  match bt with
  | SingleBond => 305.0
  | DoubleBond => 615.0
  | TripleBond => 891.0
  | _          => 0.0
  end.

Definition bond_energy_CO (bt : BondType) : R :=
  match bt with
  | SingleBond => 360.0
  | DoubleBond => 736.0
  | _          => 0.0
  end.

Definition bond_energy_CH  : R := 413.0.
Definition bond_energy_NH  : R := 391.0.
Definition bond_energy_OH  : R := 459.0.
Definition bond_energy_CF  : R := 485.0.
Definition bond_energy_CCl : R := 339.0.
Definition bond_energy_CBr : R := 285.0.
Definition bond_energy_CI  : R := 218.0.
Definition bond_energy_CS  : R := 272.0.
Definition bond_energy_NN  : R := 163.0.
Definition bond_energy_OO  : R := 146.0.

(** Bond energy increases with bond order for C-C *)
Theorem CC_energy_increases_with_order :
  bond_energy_CC SingleBond < bond_energy_CC DoubleBond /\
  bond_energy_CC DoubleBond < bond_energy_CC TripleBond.
Proof.
  unfold bond_energy_CC; simpl; lra.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 7. Bond Polarity                                                 *)
(** ------------------------------------------------------------------ *)

(** Bond polarity (difference in electronegativity between the two atoms)
    A bond is considered polar if |ΔEN| > 0.5 on the Pauling scale,
    ionic if |ΔEN| > 1.7. *)

Definition delta_electronegativity (e1 e2 : Element) : R :=
  Rabs (electronegativity e1 - electronegativity e2).

Inductive PolarityClass : Type :=
  | Nonpolar   (** |ΔEN| < 0.5   *)
  | Polar      (** 0.5 ≤ |ΔEN| ≤ 1.7 *)
  | Ionic.     (** |ΔEN| > 1.7   *)

Definition classify_polarity (e1 e2 : Element) : PolarityClass :=
  let dEN := delta_electronegativity e1 e2 in
  if Rlt_dec dEN 0.5 then Nonpolar
  else if Rlt_dec dEN 1.7 then Polar
  else Ionic.

(** ------------------------------------------------------------------ *)
(** ** 8. Steric Constraints (Van der Waals overlap)                    *)
(** ------------------------------------------------------------------ *)

(** Bonded atoms: distance ≥ 0.5 × sum of VdW radii *)
Definition min_bonded_distance (e1 e2 : Element) : R :=
  0.5 * (van_der_waals_radius e1 + van_der_waals_radius e2).

(** Non-bonded atoms in the same molecule: distance ≥ 0.8 × sum VdW *)
Definition min_nonbonded_distance (e1 e2 : Element) : R :=
  0.8 * (van_der_waals_radius e1 + van_der_waals_radius e2).

(** Atoms in different molecules (non-interacting): distance ≥ sum VdW *)
Definition min_noninteracting_distance (e1 e2 : Element) : R :=
  van_der_waals_radius e1 + van_der_waals_radius e2.

(** No atomic overlap: bonded distance threshold is less than non-bonded.
    Uses vdw_radius_positive from Atoms to avoid a 118×118 case split. *)
Theorem bonded_lt_nonbonded : forall e1 e2 : Element,
    min_bonded_distance e1 e2 < min_nonbonded_distance e1 e2.
Proof.
  intros e1 e2.
  unfold min_bonded_distance, min_nonbonded_distance.
  assert (H1 := vdw_radius_positive e1).
  assert (H2 := vdw_radius_positive e2).
  lra.
Qed.

(** Non-bonded distance threshold is less than non-interacting *)
Theorem nonbonded_lt_noninteracting : forall e1 e2 : Element,
    min_nonbonded_distance e1 e2 < min_noninteracting_distance e1 e2.
Proof.
  intros e1 e2.
  unfold min_nonbonded_distance, min_noninteracting_distance.
  assert (H1 := vdw_radius_positive e1).
  assert (H2 := vdw_radius_positive e2).
  lra.
Qed.

Close Scope R_scope.
