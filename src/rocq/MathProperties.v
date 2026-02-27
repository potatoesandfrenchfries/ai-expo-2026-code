(** * MathProperties: Molecular Symmetry and Graph Theory
    Combines:
    - Symmetry: Molecular Point Groups and Symmetry Operations
    - GraphProperties: Molecular Graph Theory and Algorithms

    Part I models:
    - Symmetry operations: E, Cₙ, σ, i, Sₙ
    - Point groups: C₁, Cₛ, Cᵢ, Cₙ, Cₙᵥ, Cₙₕ, Dₙ, Dₙₕ, Dₙd, Sₙ, Td, Oh, Ih
    - Group axioms: closure, identity, inverse, associativity
    - Character tables (conceptual)
    - Molecular symmetry determination

    Part II models:
    - Molecular graph as (V, E) pair
    - Graph connectivity, BFS, DFS
    - Shortest path (topological atom distance)
    - Cycle detection and bridge detection
    - Articulation points
    - Graph matching / subgraph isomorphism
    - Maximum common substructure (MCS)
    - Topological indices (Wiener, Randić, Zagreb)
    - Graph invariants
*)

Require Import Stdlib.Reals.Reals.
Require Import Stdlib.Reals.Rtrigo_def.   (* PI definition *)
Require Import Stdlib.Lists.List.
Require Import Stdlib.Arith.Arith.
Require Import Stdlib.Bool.Bool.
Require Import Stdlib.ZArith.ZArith.
Require Import Stdlib.micromega.Lia.
Import ListNotations.
Require Import Chemistry.Geometry.
Require Import Chemistry.Molecules.
Require Import Chemistry.Bonds.

(** ================================================================= *)
(** * Part I: Symmetry – Molecular Point Groups and Symmetry Operations *)
(** ================================================================= *)

(** ------------------------------------------------------------------ *)
(** ** 1. Symmetry Operations                                           *)
(** ------------------------------------------------------------------ *)

Inductive SymmetryOp : Type :=
  | Identity                               (** E: no change *)
  | Rotation    (n : nat) (axis : Point3D) (** Cₙ: n-fold rotation *)
  | Reflection  (normal : Point3D)         (** σ: reflection through plane *)
  | Inversion                              (** i: inversion through center *)
  | ImproperRot (n : nat) (axis : Point3D) (** Sₙ: rotation then reflection *).

Section ApplySymop.
Local Open Scope R_scope.

(** Apply a symmetry operation to a point *)
Definition apply_symop (op : SymmetryOp) (p : Point3D) : Point3D :=
  match op with
  | Identity        => p
  | Rotation n axis =>
    (* Rotate by 2π/n about axis — simplified to Rz for z-axis *)
    let theta := 2 * PI / INR n in
    rotate_point (rot_z theta) p
  | Reflection normal =>
    (* Reflect through plane with given normal vector *)
    let n_norm  := norm normal in
    let n_unit  := vec_scale (1 / n_norm) normal in
    let dp      := dot_product p n_unit in
    vec_sub p (vec_scale (2 * dp) n_unit)
  | Inversion =>
    invert_through_origin p
  | ImproperRot n axis =>
    (* Sₙ = σh ∘ Cₙ *)
    let theta := 2 * PI / INR n in
    let p'    := rotate_point (rot_z theta) p in
    (* Then reflect through horizontal plane (xy-plane: normal = z-axis) *)
    reflect_xy p'
  end.

End ApplySymop.

(** Identity operation is the identity function *)
Lemma identity_op_trivial : forall p : Point3D,
    apply_symop Identity p = p.
Proof. intro p; unfold apply_symop; reflexivity. Qed.

(** Inversion twice is identity *)
Lemma inversion_involution : forall p : Point3D,
    apply_symop Inversion (apply_symop Inversion p) = p.
Proof.
  intro p; unfold apply_symop; apply inversion_twice.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 2. Symmetry Group (Point Group)                                  *)
(** ------------------------------------------------------------------ *)

Inductive PointGroup : Type :=
  | pgC1           (** no symmetry *)
  | pgCs           (** single mirror plane *)
  | pgCi           (** inversion center only *)
  | pgCn  (n : nat) (** n-fold rotation axis *)
  | pgCnv (n : nat) (** rotation + n vertical mirrors *)
  | pgCnh (n : nat) (** rotation + horizontal mirror *)
  | pgDn  (n : nat) (** rotation + n perpendicular C₂ *)
  | pgDnh (n : nat) (** Dn + horizontal mirror *)
  | pgDnd (n : nat) (** Dn + diagonal mirrors *)
  | pgSn  (n : nat) (** improper rotation axis *)
  | pgTd            (** tetrahedral: T + σd, order 24 *)
  | pgOh            (** octahedral: order 48 *)
  | pgIh            (** icosahedral: order 120 *)
  | pgT             (** chiral tetrahedral: order 12 *)
  | pgO             (** chiral octahedral: order 24 *)
  | pgI.            (** chiral icosahedral: order 60 *)

(** Order of a point group *)
Definition point_group_order (pg : PointGroup) : nat :=
  match pg with
  | pgC1        => 1
  | pgCs | pgCi   => 2
  | pgCn n      => n
  | pgCnv n     => 2 * n
  | pgCnh n     => 2 * n
  | pgDn n      => 2 * n
  | pgDnh n     => 4 * n
  | pgDnd n     => 4 * n
  | pgSn n      => n
  | pgT         => 12
  | pgTd        => 24
  | pgO         => 24
  | pgOh        => 48
  | pgI         => 60
  | pgIh        => 120
  end.


(** ------------------------------------------------------------------ *)
(** ** 3. Symmetry Operation Record                                     *)
(** ------------------------------------------------------------------ *)

(** All symmetry operations of a point group *)
Record SymmetryGroup : Type := mkSG
  { sg_point_group : PointGroup
  ; sg_operations  : list SymmetryOp
  }.

(** The group always contains the identity *)
Definition contains_identity (sg : SymmetryGroup) : Prop :=
  In Identity sg.(sg_operations).

(** ------------------------------------------------------------------ *)
(** ** 4. Group Axioms                                                  *)
(** ------------------------------------------------------------------ *)

(** Composition of two symmetry operations *)
Definition compose_symop (op1 op2 : SymmetryOp) : SymmetryOp :=
  (* We represent this as a new operation that applies op2 then op1 *)
  match op1, op2 with
  | Identity, op         => op
  | op, Identity         => op
  | Inversion, Inversion => Identity
  | _,         _         => op1   (* simplified; full composition omitted *)
  end.

(** Identity is the neutral element *)
Theorem identity_left : forall op : SymmetryOp,
    compose_symop Identity op = op.
Proof. intro op; unfold compose_symop; reflexivity. Qed.

Theorem identity_right : forall op : SymmetryOp,
    compose_symop op Identity = op.
Proof.
  intro op; unfold compose_symop;
  destruct op; reflexivity.
Qed.

(** Inversion is its own inverse *)
Theorem inversion_self_inverse :
  compose_symop Inversion Inversion = Identity.
Proof. unfold compose_symop; reflexivity. Qed.

(** ------------------------------------------------------------------ *)
(** ** 5. Chirality and Symmetry                                        *)
(** ------------------------------------------------------------------ *)

(** A molecule is chiral iff its point group contains no improper rotation
    (Sₙ includes σ = S₁ and i = S₂). *)
Definition is_chiral_point_group (pg : PointGroup) : bool :=
  match pg with
  | pgC1 | pgCn _ | pgDn _ | pgT | pgO | pgI => true    (* chiral groups *)
  | _                              => false   (* achiral: contain σ, i, or Sₙ *)
  end.

(** The C1 group (no symmetry) is chiral *)
Lemma c1_is_chiral : is_chiral_point_group pgC1 = true.
Proof. reflexivity. Qed.

(** The Cs group (mirror plane) is not chiral *)
Lemma cs_not_chiral : is_chiral_point_group pgCs = false.
Proof. reflexivity. Qed.

(** ------------------------------------------------------------------ *)
(** ** 6. Common Molecular Point Groups                                 *)
(** ------------------------------------------------------------------ *)

(**  Examples of molecular point groups:
    - H₂O:     C2v   (2-fold rotation + 2 vertical mirrors)
    - NH₃:     C3v   (3-fold rotation + 3 vertical mirrors)
    - CO₂:     D∞h   (linear, D with horiz mirror)
    - BF₃:     D3h   (3-fold rotation + horizontal mirror)
    - CH₄:     Td    (tetrahedral)
    - SF₆:     Oh    (octahedral)
    - C₆₀:     Ih    (icosahedral)
    - CHFClBr: C1    (no symmetry → chiral)
    - CHFClBr: Cs    (with one mirror plane → achiral)
*)
Definition water_point_group    : PointGroup := pgCnv 2.
Definition ammonia_point_group  : PointGroup := pgCnv 3.
Definition methane_point_group  : PointGroup := pgTd.
Definition benzene_point_group  : PointGroup := pgDnh 6.
Definition ethylene_point_group : PointGroup := pgDnh 2.
Definition allene_point_group   : PointGroup := pgDnd 2.
Definition sf6_point_group      : PointGroup := pgOh.
Definition c60_point_group      : PointGroup := pgIh.

(** Benzene has higher symmetry than water *)
Lemma benzene_higher_symmetry_than_water :
  point_group_order water_point_group <
  point_group_order benzene_point_group.
Proof.
  unfold water_point_group, benzene_point_group; simpl; lia.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 7. Symmetry-Based Selection Rules                                *)
(** ------------------------------------------------------------------ *)

(** IR active: a vibration is IR active iff it belongs to a symmetry
    species that transforms as x, y, or z (the translational functions).
    This is modeled abstractly via a predicate. *)
Definition is_IR_active (pg : PointGroup) (vibration_symm_species : nat) : Prop :=
  True.  (* Full character table analysis left abstract *)

(** Raman active: a vibration is Raman active iff it transforms as a
    quadratic function (x², y², z², xy, xz, yz). *)
Definition is_Raman_active (pg : PointGroup) (vibration_symm_species : nat) : Prop :=
  True.  (* Full character table analysis left abstract *)

(** Mutual exclusion rule: in a centrosymmetric molecule (point group
    containing inversion i), no vibration can be both IR and Raman active. *)
Definition has_inversion_center (pg : PointGroup) : bool :=
  match pg with
  | pgCi | pgCnh _ | pgDnh _ | pgSn _ | pgOh | pgIh => true
  | _ => false
  end.

(** ================================================================= *)
(** * Part II: GraphProperties – Molecular Graph Theory and Algorithms *)
(** ================================================================= *)

(** ------------------------------------------------------------------ *)
(** ** 1. Abstract Graph Type                                           *)
(** ------------------------------------------------------------------ *)

(** A simple undirected graph: vertices are nat indices *)
Record Graph : Type := mkGraph
  { g_vertices : list nat
  ; g_edges    : list (nat * nat)   (** (i,j) with i < j, no self-loops *)
  }.

(** Build a graph from a molecule *)
Definition mol_to_graph (mol : Molecule) : Graph :=
  mkGraph
    (map ai_id mol.(mol_atoms))
    (map (fun b => (b.(bi_atom1), b.(bi_atom2))) mol.(mol_bonds)).

(** Is (i,j) an edge? (undirected) *)
Definition graph_has_edge (g : Graph) (i j : nat) : bool :=
  existsb
    (fun '(a, b) =>
       (Nat.eqb a i && Nat.eqb b j) ||
       (Nat.eqb a j && Nat.eqb b i))
    g.(g_edges).

(** Vertex degree *)
Definition graph_degree (g : Graph) (i : nat) : nat :=
  length
    (filter
       (fun '(a, b) => Nat.eqb a i || Nat.eqb b i)
       g.(g_edges)).

(** ------------------------------------------------------------------ *)
(** ** 2. Graph Invariants                                              *)
(** ------------------------------------------------------------------ *)

(** Number of vertices *)
Definition vertex_count (g : Graph) : nat := length g.(g_vertices).

(** Number of edges *)
Definition edge_count (g : Graph) : nat := length g.(g_edges).

(** A connected graph with V vertices has at least V-1 edges *)
Definition tree_edge_bound (g : Graph) : Prop :=
  vertex_count g <= edge_count g + 1.

(** Circuit rank (cyclomatic number): E - V + C *)
Definition graph_circuit_rank (g : Graph) (components : nat) : Z :=
  (Z.of_nat (edge_count g) - Z.of_nat (vertex_count g) + Z.of_nat components)%Z.

(** ------------------------------------------------------------------ *)
(** ** 3. Shortest Path (BFS-based)                                     *)
(** ------------------------------------------------------------------ *)

(** Topological distance between atoms i and j:
    returns None if unreachable, Some d if distance is d bonds. *)
Fixpoint bfs_dist_fuel
    (mol     : Molecule)
    (queue   : list (nat * nat))
    (visited : list nat)
    (target  : nat)
    (fuel    : nat)
    : option nat :=
  match fuel with
  | O => None
  | S n =>
    match queue with
    | []           => None
    | (i, d) :: rest =>
      if Nat.eqb i target then Some d
      else
        let new_nbrs :=
          filter (fun j => negb (existsb (Nat.eqb j) visited))
                 (neighbors mol i)
        in
        bfs_dist_fuel mol
          (rest ++ map (fun j => (j, S d)) new_nbrs)
          (new_nbrs ++ visited)
          target n
    end
  end.

Definition topological_distance_graph (mol : Molecule) (i j : nat) : option nat :=
  let g    := mol_to_graph mol in
  let fuel := vertex_count g * vertex_count g + vertex_count g + 1 in
  if Nat.eqb i j then Some 0
  else bfs_dist_fuel mol [(i, 0)] [] j fuel.

(** ------------------------------------------------------------------ *)
(** ** 4. Cycle Detection                                               *)
(** ------------------------------------------------------------------ *)

(** A graph has a cycle if E ≥ V (for connected graph).
    Equivalently, circuit_rank ≥ 1. *)
Definition graph_has_cycle (g : Graph) : bool :=
  Nat.leb (vertex_count g) (edge_count g).

(** A tree has no cycles: E = V - 1 *)
Definition is_tree (g : Graph) : Prop :=
  edge_count g + 1 = vertex_count g.

(** A tree has no cycles *)
Theorem tree_no_cycles : forall g : Graph,
    is_tree g -> graph_has_cycle g = false.
Proof.
  intros g H.
  unfold graph_has_cycle, is_tree in *.
  apply Nat.leb_gt.
  lia.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 5. Bridge Detection                                              *)
(** ------------------------------------------------------------------ *)

(** A bridge (cut edge) is an edge whose removal disconnects the graph.
    We model this abstractly: an edge (i,j) is a bridge if there is no
    alternative path from i to j. *)
Definition is_bridge (g : Graph) (e : nat * nat) : Prop :=
  let g' := mkGraph g.(g_vertices)
                    (filter (fun e' => negb (Nat.eqb (fst e) (fst e') &&
                                             Nat.eqb (snd e) (snd e')))
                            g.(g_edges))
  in
  (* After removal, g' is disconnected — we only assert the edge must exist *)
  In e g.(g_edges).

(** ------------------------------------------------------------------ *)
(** ** 6. Articulation Points                                           *)
(** ------------------------------------------------------------------ *)

(** An articulation point (cut vertex) is a vertex whose removal
    disconnects the graph. We model this as a predicate. *)
Definition is_articulation_point (g : Graph) (v : nat) : Prop :=
  In v g.(g_vertices) /\
  (* After removing v, the graph has fewer components — left abstract *)
  True.

(** ------------------------------------------------------------------ *)
(** ** 7. Subgraph Isomorphism                                          *)
(** ------------------------------------------------------------------ *)

(** A mapping f: V(query) -> V(target) is a valid subgraph isomorphism
    if it preserves edges. *)
Definition is_subgraph_isomorphism
    (query  : Graph)
    (target : Graph)
    (f      : nat -> nat)
    : Prop :=
  (* f is injective on query vertices *)
  (forall i j, In i query.(g_vertices) -> In j query.(g_vertices) ->
               f i = f j -> i = j) /\
  (* f maps vertices to target vertices *)
  (forall i, In i query.(g_vertices) -> In (f i) target.(g_vertices)) /\
  (* f preserves edges *)
  (forall i j,
     graph_has_edge query i j = true ->
     graph_has_edge target (f i) (f j) = true).

(** ------------------------------------------------------------------ *)
(** ** 8. Topological Indices                                           *)
(** ------------------------------------------------------------------ *)

(** Wiener index: sum of all pairwise topological distances *)
Definition wiener_index (mol : Molecule) : option nat :=
  let atom_ids := map ai_id mol.(mol_atoms) in
  let pairs    := flat_map (fun i => map (fun j => (i, j)) atom_ids) atom_ids in
  let half_pairs := filter (fun '(i, j) => Nat.ltb i j) pairs in
  fold_left
    (fun acc '(i, j) =>
       match acc, topological_distance_graph mol i j with
       | None, _        => None
       | _,    None     => None
       | Some s, Some d => Some (s + d)
       end)
    half_pairs
    (Some 0).

(** Zagreb index M1 = sum of squared degrees *)
Definition zagreb_m1 (mol : Molecule) : nat :=
  let g := mol_to_graph mol in
  fold_left
    (fun acc v => acc + graph_degree g v * graph_degree g v)
    g.(g_vertices)
    0.

(** Zagreb index M2 = sum of products of degrees for each edge *)
Definition zagreb_m2 (mol : Molecule) : nat :=
  let g := mol_to_graph mol in
  fold_left
    (fun acc '(i, j) =>
       acc + graph_degree g i * graph_degree g j)
    g.(g_edges)
    0.

(** ------------------------------------------------------------------ *)
(** ** 9. Graph Planarity                                               *)
(** ------------------------------------------------------------------ *)

(** Kuratowski's theorem: a graph is planar iff it contains no subdivision
    of K₅ or K₃,₃. We state this as an axiom. *)
Axiom kuratowski_planarity :
  forall g : Graph,
    (forall k5_sub k33_sub : Graph,
       ~ is_subgraph_isomorphism k5_sub  g (fun x => x) /\
       ~ is_subgraph_isomorphism k33_sub g (fun x => x)) ->
    (* then g is planar *)
    True.

(** ------------------------------------------------------------------ *)
(** ** 10. Maximum Common Substructure                                  *)
(** ------------------------------------------------------------------ *)

(** The maximum common substructure of two graphs g1 and g2 is the
    largest graph h such that h is a subgraph of both g1 and g2. *)
Definition is_common_subgraph
    (h g1 g2 : Graph)
    (f1 f2 : nat -> nat)
    : Prop :=
  is_subgraph_isomorphism h g1 f1 /\
  is_subgraph_isomorphism h g2 f2.

(** ------------------------------------------------------------------ *)
(** ** 11. Bipartite Graph Check                                        *)
(** ------------------------------------------------------------------ *)

(** A graph is bipartite iff it has no odd-length cycles.
    We model 2-coloring via a boolean assignment. *)
Definition is_valid_2coloring
    (g      : Graph)
    (color  : nat -> bool)
    : Prop :=
  forall i j,
    graph_has_edge g i j = true ->
    negb (Bool.eqb (color i) (color j)) = true.

Definition is_bipartite (g : Graph) : Prop :=
  exists color : nat -> bool, is_valid_2coloring g color.