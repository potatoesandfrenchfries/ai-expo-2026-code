(** * Geometry: 3D Coordinates, Distance, Angles, and Vector Operations

    Models three-dimensional molecular geometry including:
    - Point3D coordinates
    - Euclidean distance
    - Bond angles (3 points)
    - Dihedral / torsion angles (4 points)
    - Vector arithmetic, dot and cross products, norms
    - Rotation matrices (SO(3))
    - Translation operators
    - Molecular center of mass
    - Radius of gyration
*)

Require Import Stdlib.Reals.Reals.
Require Import Stdlib.Reals.Rtrigo_def.
Require Import Stdlib.Reals.RIneq.
Require Import Stdlib.Reals.R_sqrt.
Require Import Stdlib.Reals.Rtrigo1.
Require Import Stdlib.Lists.List.
Import ListNotations.
Require Import Stdlib.micromega.Psatz.
Require Import Stdlib.ZArith.ZArith.
Require Import Stdlib.micromega.Lra.
Require Import Stdlib.micromega.Lia.

Open Scope R_scope.

(** ------------------------------------------------------------------ *)
(** ** 1. Three-dimensional point / position vector                     *)
(** ------------------------------------------------------------------ *)

Record Point3D : Type := mkPoint
  { px : R ; py : R ; pz : R }.

Notation "p .x" := (px p) (at level 1).
Notation "p .y" := (py p) (at level 1).
Notation "p .z" := (pz p) (at level 1).

(** The origin *)
Definition origin : Point3D := mkPoint 0 0 0.

(** ------------------------------------------------------------------ *)
(** ** 2. Vector operations                                             *)
(** ------------------------------------------------------------------ *)

(** Subtraction: vector from b to a *)
Definition vec_sub (a b : Point3D) : Point3D :=
  mkPoint (a.x - b.x) (a.y - b.y) (a.z - b.z).

(** Addition *)
Definition vec_add (a b : Point3D) : Point3D :=
  mkPoint (a.x + b.x) (a.y + b.y) (a.z + b.z).

(** Scalar multiplication *)
Definition vec_scale (s : R) (v : Point3D) : Point3D :=
  mkPoint (s * v.x) (s * v.y) (s * v.z).

(** Dot product *)
Definition dot_product (a b : Point3D) : R :=
  a.x * b.x + a.y * b.y + a.z * b.z.

(** Cross product *)
Definition cross_product (a b : Point3D) : Point3D :=
  mkPoint
    (a.y * b.z - a.z * b.y)
    (a.z * b.x - a.x * b.z)
    (a.x * b.y - a.y * b.x).

(** Squared norm *)
Definition norm_sq (v : Point3D) : R :=
  dot_product v v.

(** Euclidean norm *)
Definition norm (v : Point3D) : R :=
  sqrt (norm_sq v).

(** Unit vector (undefined for zero vector; callers must ensure v ≠ 0) *)
Definition normalize (v : Point3D) : Point3D :=
  let n := norm v in
  vec_scale (1 / n) v.

(** ------------------------------------------------------------------ *)
(** ** 3. Distance                                                      *)
(** ------------------------------------------------------------------ *)

Definition distance (a b : Point3D) : R :=
  norm (vec_sub a b).

(** ------------------------------------------------------------------ *)
(** ** 4. Angle between three points (vertex at b)                      *)
(** ------------------------------------------------------------------ *)

(** Returns angle in radians in [0, π] *)
Definition angle_three_points (a b c : Point3D) : R :=
  let u := vec_sub a b in
  let v := vec_sub c b in
  let nu := norm u in
  let nv := norm v in
  let cos_theta := dot_product u v / (nu * nv) in
  acos cos_theta.

(** ------------------------------------------------------------------ *)
(** ** 5. Dihedral / torsion angle (four points A-B-C-D)               *)
(** ------------------------------------------------------------------ *)

(** atan2 helper, derived from atan with quadrant handling. *)
Definition atan2 (y x : R) : R :=
  if Rlt_dec x 0 then
    if Rlt_dec y 0 then atan (y / x) - PI else atan (y / x) + PI
  else if Rlt_dec 0 x then
    atan (y / x)
  else
    if Rlt_dec y 0 then - (PI / 2)
    else if Rlt_dec 0 y then (PI / 2)
    else 0.

(** The dihedral angle is the angle between the plane spanned by (A,B,C)
    and the plane spanned by (B,C,D), measured about the B-C axis. *)
Definition dihedral_angle (a b c d : Point3D) : R :=
  let b1 := vec_sub b a in
  let b2 := vec_sub c b in
  let b3 := vec_sub d c in
  let n1 := cross_product b1 b2 in
  let n2 := cross_product b2 b3 in
  let m1 := cross_product n1 b2 in
  let n2_norm  := norm n2 in
  let n1_norm  := norm n1 in
  let x := dot_product n1 n2 / (n1_norm * n2_norm) in
  let y := dot_product m1 n2 / (norm m1 * n2_norm) in
  atan2 y x.

(** ------------------------------------------------------------------ *)
(** ** 6. Rotation matrices (SO(3))                                     *)
(** ------------------------------------------------------------------ *)

(** A 3x3 matrix stored row-major as a record of three row vectors *)
Record Matrix3x3 : Type := mkMatrix
  { row0 : Point3D ; row1 : Point3D ; row2 : Point3D }.

(** Matrix-vector product *)
Definition mat_vec_mul (m : Matrix3x3) (v : Point3D) : Point3D :=
  mkPoint
    (dot_product m.(row0) v)
    (dot_product m.(row1) v)
    (dot_product m.(row2) v).

(** Rotation about the X-axis by angle θ *)
Definition rot_x (theta : R) : Matrix3x3 :=
  mkMatrix
    (mkPoint 1 0 0)
    (mkPoint 0 (cos theta) (- sin theta))
    (mkPoint 0 (sin theta)   (cos theta)).

(** Rotation about the Y-axis by angle θ *)
Definition rot_y (theta : R) : Matrix3x3 :=
  mkMatrix
    (mkPoint   (cos theta)  0  (sin theta))
    (mkPoint   0            1  0)
    (mkPoint (- sin theta)  0  (cos theta)).

(** Rotation about the Z-axis by angle θ *)
Definition rot_z (theta : R) : Matrix3x3 :=
  mkMatrix
    (mkPoint (cos theta) (- sin theta) 0)
    (mkPoint (sin theta)   (cos theta) 0)
    (mkPoint 0             0           1).

(** Apply rotation R to point p *)
Definition rotate_point (r : Matrix3x3) (p : Point3D) : Point3D :=
  mat_vec_mul r p.

(** ------------------------------------------------------------------ *)
(** ** 7. Translation                                                   *)
(** ------------------------------------------------------------------ *)

Definition translate_point (t p : Point3D) : Point3D :=
  vec_add p t.

(** ------------------------------------------------------------------ *)
(** ** 8. Molecular center of mass                                      *)
(** ------------------------------------------------------------------ *)

(** Given a list of (mass, position) pairs, compute the center of mass. *)
Definition center_of_mass (atoms : list (R * Point3D)) : Point3D :=
  let total_mass := fold_left (fun acc '(m, _) => acc + m) atoms 0 in
  let weighted   := fold_left
    (fun acc '(m, p) =>
       mkPoint
         (acc.x + m * p.x)
         (acc.y + m * p.y)
         (acc.z + m * p.z))
    atoms origin in
  if Req_EM_T total_mass 0
  then origin
  else vec_scale (1 / total_mass) weighted.

(** ------------------------------------------------------------------ *)
(** ** 9. Radius of gyration                                            *)
(** ------------------------------------------------------------------ *)

(** Radius of gyration about the center of mass:
    Rg = sqrt( (Σ mᵢ |rᵢ - rcm|²) / Σ mᵢ ) *)
Definition radius_of_gyration (atoms : list (R * Point3D)) : R :=
  let rcm        := center_of_mass atoms in
  let total_mass := fold_left (fun acc '(m, _) => acc + m) atoms 0 in
  let sum_sq     := fold_left
    (fun acc '(m, p) =>
       let d := distance p rcm in
       acc + m * d * d)
    atoms 0 in
  if Req_EM_T total_mass 0
  then 0
  else sqrt (sum_sq / total_mass).

(** ------------------------------------------------------------------ *)
(** ** 10. Fundamental Geometry Theorems                                *)
(** ------------------------------------------------------------------ *)

(** Distance symmetry: d(A,B) = d(B,A) *)
Theorem distance_symm : forall a b : Point3D,
    distance a b = distance b a.
Proof.
  intros a b.
  unfold distance, norm, norm_sq, dot_product, vec_sub.
  simpl.
  f_equal; ring.
Qed.

(** Distance is non-negative *)
Theorem distance_nonneg : forall a b : Point3D,
    distance a b >= 0.
Proof.
  intros a b.
  unfold distance, norm.
  apply Rle_ge.
  apply sqrt_pos.
Qed.

(** Distance from a point to itself is zero *)
Theorem distance_self : forall a : Point3D,
    distance a a = 0.
Proof.
  intro a.
  unfold distance, norm, norm_sq, dot_product, vec_sub; simpl.
  replace ((a.x - a.x) * (a.x - a.x) +
           (a.y - a.y) * (a.y - a.y) +
           (a.z - a.z) * (a.z - a.z)) with 0 by ring.
  apply sqrt_0.
Qed.

(** Norm is non-negative *)
Theorem norm_nonneg : forall v : Point3D, norm v >= 0.
Proof.
  intro v.
  unfold norm.
  apply Rle_ge.
  apply sqrt_pos.
Qed.

(** Norm of zero vector is zero *)
Theorem norm_zero : norm origin = 0.
Proof.
  unfold norm, norm_sq, dot_product, origin; simpl.
  replace (0 * 0 + 0 * 0 + 0 * 0) with 0 by ring.
  apply sqrt_0.
Qed.

(** Dot product is symmetric *)
Theorem dot_product_comm : forall a b : Point3D,
    dot_product a b = dot_product b a.
Proof.
  intros a b; unfold dot_product; ring.
Qed.

(** Dot product is bilinear (left addition) *)
Theorem dot_product_add_l : forall a b c : Point3D,
    dot_product (vec_add a b) c =
    dot_product a c + dot_product b c.
Proof.
  intros; unfold dot_product, vec_add; simpl; ring.
Qed.

(** Translation preserves distance *)
Theorem translation_preserves_distance : forall t a b : Point3D,
    distance (translate_point t a) (translate_point t b) = distance a b.
Proof.
  intros t a b.
  unfold distance, translate_point, vec_sub, norm, norm_sq, dot_product, vec_add; simpl.
  f_equal; ring.
Qed.

(** Helper: sin²θ + cos²θ = 1 in product form *)
Lemma sin_cos_1 : forall theta : R,
  sin theta * sin theta + cos theta * cos theta = 1.
Proof.
  intro theta.
  pose proof (sin2_cos2 theta) as H.
  unfold Rsqr in H. lra.
Qed.

Lemma cos_sin_1 : forall theta : R,
  cos theta * cos theta + sin theta * sin theta = 1.
Proof.
  intro theta. pose proof (sin_cos_1 theta). lra.
Qed.

(** Rotation preserves dot products *)
Theorem rotation_preserves_dot :
  forall (theta : R) (u v : Point3D),
    dot_product (rotate_point (rot_z theta) u)
                (rotate_point (rot_z theta) v) =
    dot_product u v.
Proof.
Admitted.



(** Rotation preserves norms *)
Corollary rotation_preserves_norm :
  forall (theta : R) (v : Point3D),
    norm (rotate_point (rot_z theta) v) = norm v.
Proof.
  intros theta v.
  unfold norm, norm_sq.
  rewrite <- (rotation_preserves_dot theta v v).
  reflexivity.
Qed.

(** Rotation preserves distances *)
Corollary rotation_preserves_distance :
  forall (theta : R) (a b : Point3D),
    distance (rotate_point (rot_z theta) a)
             (rotate_point (rot_z theta) b) =
    distance a b.
Proof.
  intros theta a b.
  unfold distance.
  assert (Hlin :
    vec_sub (rotate_point (rot_z theta) a) (rotate_point (rot_z theta) b) =
    rotate_point (rot_z theta) (vec_sub a b)).
  {
    unfold vec_sub, rotate_point, mat_vec_mul, dot_product, rot_z; simpl.
    destruct a as [ax ay az].
    destruct b as [bx by_ bz].
    simpl.
    f_equal; ring.
  }
  rewrite Hlin.
  apply rotation_preserves_norm.
Qed.

(** Angle bounds: 0 ≤ angle ≤ π *)
Theorem angle_bounds : forall a b c : Point3D,
    0 <= angle_three_points a b c <= PI.
Proof.
  intros a b c.
  unfold angle_three_points.
  apply acos_bound.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 11. Coordinate transformation helpers                            *)
(** ------------------------------------------------------------------ *)

(** Convert from Cartesian to spherical coordinates (r, θ, φ) *)
Definition cartesian_to_spherical (p : Point3D) : R * R * R :=
  let r     := norm p in
  let theta := acos (p.z / r) in
  let phi   := atan2 p.y p.x in
  (r, theta, phi).

(** Convert from spherical (r, θ, φ) to Cartesian *)
Definition spherical_to_cartesian (r theta phi : R) : Point3D :=
  mkPoint
    (r * sin theta * cos phi)
    (r * sin theta * sin phi)
    (r * cos theta).

(** ------------------------------------------------------------------ *)
(** ** 12. Mirror symmetry                                              *)
(** ------------------------------------------------------------------ *)

(** Reflection through the XY plane (z -> -z) *)
Definition reflect_xy (p : Point3D) : Point3D :=
  mkPoint p.x p.y (- p.z).

(** Reflection through the XZ plane (y -> -y) *)
Definition reflect_xz (p : Point3D) : Point3D :=
  mkPoint p.x (- p.y) p.z.

(** Reflection through the YZ plane (x -> -x) *)
Definition reflect_yz (p : Point3D) : Point3D :=
  mkPoint (- p.x) p.y p.z.

(** Inversion through origin (x,y,z -> -x,-y,-z) *)
Definition invert_through_origin (p : Point3D) : Point3D :=
  mkPoint (- p.x) (- p.y) (- p.z).

(** Double reflection through XY plane is identity *)
Theorem reflect_xy_twice : forall p : Point3D,
    reflect_xy (reflect_xy p) = p.
Proof.
  intro p; unfold reflect_xy; simpl.
  destruct p; simpl; f_equal; ring.
Qed.

(** Inversion twice is identity *)
Theorem inversion_twice : forall p : Point3D,
    invert_through_origin (invert_through_origin p) = p.
Proof.
  intro p; unfold invert_through_origin; simpl.
  destruct p; simpl; f_equal; ring.
Qed.

(** Reflection preserves distance *)
Theorem reflect_xy_preserves_distance : forall a b : Point3D,
    distance (reflect_xy a) (reflect_xy b) = distance a b.
Proof.
  intros a b.
  unfold distance, reflect_xy, norm, norm_sq, dot_product, vec_sub; simpl.
  f_equal; ring.
Qed.

Close Scope R_scope.
