(** * HydrogenBonding: H-Bond Criteria, Types, and Constraints
    Models:
    - H-bond donor elements (N-H, O-H, F-H)
    - H-bond acceptor elements (N, O, F with lone pairs)
    - Distance criterion: heavy-atom distance 2.5-3.5 Å
    - Angle criterion: D-H...A ≥ 120°
    - H-bond strength: 1-5 kcal/mol
    - H-bond types: intra/intermolecular, bifurcated, cooperative
*)
 
Require Import Stdlib.Reals.Reals.
Require Import Stdlib.Reals.Rtrigo_def.   (* PI *)
Require Import Stdlib.Reals.Rtrigo1.
Require Import Stdlib.Lists.List.
Require Import Stdlib.micromega.Lia.
Require Import Stdlib.micromega.Lra.
Import ListNotations.
Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.
 
Open Scope R_scope.
 
(** ------------------------------------------------------------------ *)
(** ** 1. H-Bond Donor Elements                                         *)
(** ------------------------------------------------------------------ *)
 
(** An element can donate an H-bond when bonded to hydrogen *)
Definition is_hbond_donor_element (e : Element) : bool :=
  match e with
  | eN | eO | eF => true
  | eS         => true   (* weaker donor *)
  | _         => false
  end.
 
(** ------------------------------------------------------------------ *)
(** ** 2. H-Bond Acceptor Elements                                      *)
(** ------------------------------------------------------------------ *)
 
(** An element can accept an H-bond if it has lone pairs *)
Definition is_hbond_acceptor_element (e : Element) : bool :=
  match e with
  | eN | eO | eF => true
  | eS | eCl    => true   (* weaker acceptors *)
  | _         => false
  end.
 
(** ------------------------------------------------------------------ *)
(** ** 3. H-Bond Distance Constraints (heavy-atom distance, Å)          *)
(** ------------------------------------------------------------------ *)
 
Definition hbond_dist_min : R := 2.5.
Definition hbond_dist_max : R := 3.5.
 
(** D-H distance (covalent bond to hydrogen) *)
Definition DH_dist_min : R := 0.8.
Definition DH_dist_max : R := 1.2.
 
(** H...A distance (hydrogen to acceptor) *)
Definition HA_dist_min : R := 1.5.
Definition HA_dist_max : R := 2.5.
 
(** H-bond distance constraint: heavy-atom separation is in range *)
Definition valid_hbond_distance (d_da : R) : Prop :=
  hbond_dist_min <= d_da /\ d_da <= hbond_dist_max.
 
(** ------------------------------------------------------------------ *)
(** ** 4. H-Bond Angle Constraint                                       *)
(** ------------------------------------------------------------------ *)
 
(** The D-H...A angle should be ≥ 120° (more linear = stronger H-bond) *)
Definition hbond_angle_min_deg : R := 120.0.
Definition hbond_angle_min_rad : R := 120.0 * PI / 180.
 
Definition valid_hbond_angle (theta_rad : R) : Prop :=
  theta_rad >= hbond_angle_min_rad.
 
(** A more linear H-bond (≥ 150°) is a "strong" H-bond by geometry *)
Definition strong_hbond_angle_min_rad : R := 150.0 * PI / 180.
 
Definition strong_hbond_by_angle (theta_rad : R) : Prop :=
  theta_rad >= strong_hbond_angle_min_rad.
 
(** Strong angle implies valid angle.
  Requires PI_RGT_0 so that lra can determine the sign of PI when
  comparing 150 * PI / 180 ≥ 120 * PI / 180. *)
Theorem strong_implies_valid_angle : forall theta_rad : R,
    strong_hbond_by_angle theta_rad -> valid_hbond_angle theta_rad.
Proof.
  intros theta H.
  unfold strong_hbond_by_angle in H.
  unfold valid_hbond_angle.
  unfold hbond_angle_min_rad, strong_hbond_angle_min_rad in *.
  pose proof PI_RGT_0; lra.
Qed.
 
(** ------------------------------------------------------------------ *)
(** ** 5. H-Bond Record                                                 *)
(** ------------------------------------------------------------------ *)
 
Record HBond : Type := mkHBond
  { hb_donor_atom    : nat   (* atom ID of D in D-H...A *)
  ; hb_hydrogen_atom : nat   (* atom ID of H *)
  ; hb_acceptor_atom : nat   (* atom ID of A *)
  ; hb_DA_distance   : R     (* D...A distance in Å *)
  ; hb_DHA_angle     : R     (* D-H...A angle in radians *)
  }.
 
(** A geometrically valid H-bond satisfies both distance and angle criteria *)
Definition valid_hbond (hb : HBond) : Prop :=
  valid_hbond_distance hb.(hb_DA_distance) /\
  valid_hbond_angle    hb.(hb_DHA_angle).
 
(** ------------------------------------------------------------------ *)
(** ** 6. H-Bond Strength (kcal/mol)                                    *)
(** ------------------------------------------------------------------ *)
 
Definition hbond_energy_min : R := 1.0.
Definition hbond_energy_max : R := 5.0.
 
(** Approximate H-bond energy from geometry (simplified linear model).
    Energy scales linearly from max at d = dist_min down to 0 at dist_max.
    Outside the valid range the extremes are clamped. *)
Definition estimated_hbond_energy (dist_ang : R) : R :=
  if Rlt_dec dist_ang hbond_dist_min then hbond_energy_max
  else if Rlt_dec hbond_dist_max dist_ang then 0
  else
    hbond_energy_max * (hbond_dist_max - dist_ang) /
                       (hbond_dist_max - hbond_dist_min).
 
(** ------------------------------------------------------------------ *)
(** ** 7. H-Bond Types                                                  *)
(** ------------------------------------------------------------------ *)
 
Inductive HBondType : Type :=
  | Intramolecular    (* D and A on the same molecule *)
  | Intermolecular    (* D and A on different molecules *)
  | WaterBridge       (* mediated by a bridging water molecule *)
  | Bifurcated        (* one donor hydrogen to two acceptors *)
  | Cooperative.      (* chains of H-bonds reinforcing each other *)
 
(** ------------------------------------------------------------------ *)
(** ** 8. Intramolecular H-Bond Constraints                             *)
(** ------------------------------------------------------------------ *)
 
(** Intramolecular H-bonds require the D and A to be separated by
    at least 4 bonds (1,4 interaction minimum) to form a ring *)
Definition intramolecular_hbond_min_separation : nat := 4.
 
(** ------------------------------------------------------------------ *)
(** ** 9. Bifurcated H-Bonds                                            *)
(** ------------------------------------------------------------------ *)
 
(** A bifurcated H-bond has one D-H and two acceptors A1, A2 *)
Record BifurcatedHBond : Type := mkBHB
  { bhb_donor_atom     : nat
  ; bhb_hydrogen_atom  : nat
  ; bhb_acceptor1      : nat
  ; bhb_acceptor2      : nat
  ; bhb_DA1_distance   : R
  ; bhb_DA2_distance   : R
  ; bhb_DHA1_angle     : R
  ; bhb_DHA2_angle     : R
  }.
 
Definition valid_bifurcated_hbond (bhb : BifurcatedHBond) : Prop :=
  valid_hbond_distance bhb.(bhb_DA1_distance) /\
  valid_hbond_distance bhb.(bhb_DA2_distance) /\
  valid_hbond_angle    bhb.(bhb_DHA1_angle)   /\
  valid_hbond_angle    bhb.(bhb_DHA2_angle).
 
(** ------------------------------------------------------------------ *)
(** ** 10. H-Bond Properties                                            *)
(** ------------------------------------------------------------------ *)
 
(** H-bond distance range is well-formed *)
Lemma hbond_dist_range_wf : hbond_dist_min < hbond_dist_max.
Proof. unfold hbond_dist_min, hbond_dist_max; lra. Qed.
 
(** H-bond angle threshold lies strictly between 0 and π.
  Both goals are linear in PI once PI_RGT_0 is in context. *)
Lemma hbond_angle_in_range : 0 < hbond_angle_min_rad < PI.
Proof.
  unfold hbond_angle_min_rad.
  pose proof PI_RGT_0; lra.
Qed.
 
Close Scope R_scope.