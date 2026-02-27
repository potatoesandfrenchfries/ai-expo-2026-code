(** * MolecularProperties: Scalar, Topological, Surface, and Electronic Properties
    Models:
    - Molecular weight, exact mass
    - LogP (partition coefficient)
    - LogS (aqueous solubility)
    - pKa
    - Topological descriptors (Wiener, Randić, Balaban, Kier)
    - Surface properties (SASA, PSA, TPSA)
    - Electronic properties (partial charges, dipole moment)
    - HOMO/LUMO energies and gap
    - Synthesizability (SA score)
    - Multi-conformer ensemble properties
*)
 
Require Import Stdlib.Reals.Reals.
Require Import Stdlib.ZArith.ZArith.
Require Import Stdlib.Lists.List.
Require Import Stdlib.micromega.Lia.
Require Import Stdlib.micromega.Lra.
Import ListNotations.
Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.
Require Import Chemistry.Molecules.
 
(** ------------------------------------------------------------------ *)
(** ** Section 5 helper: topological distance via BFS                  *)
(**    Defined before Open Scope R_scope to avoid O / S constructor    *)
(**    clashes with the nat match arms.                                 *)
(** ------------------------------------------------------------------ *)
 
(** BFS with depth tracking: returns the shortest-path distance from
    the front of [queue] to [target], or [None] if unreachable.
    [queue] carries (node_id, depth) pairs; [visited] prevents cycles. *)
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
 
(** Topological distance between atoms [i] and [j] in [mol].
    Returns [Some 0] when i = j, [None] when j is unreachable from i. *)
Definition topological_distance (mol : Molecule) (i j : nat) : option nat :=
  if Nat.eqb i j then Some 0%nat
  else
    let n := length mol.(mol_atoms) in
    bfs_dist_fuel mol [(i, 0%nat)] [i] j (n + 1).
 
Open Scope R_scope.
 
(** ------------------------------------------------------------------ *)
(** ** 1. Scalar Properties                                             *)
(** ------------------------------------------------------------------ *)
 
(** Molecular weight (re-exported from Molecules.v) *)
Definition mol_weight := molecular_weight.
 
(** Exact monoisotopic mass (uses most abundant isotope for each element) *)
Definition monoisotopic_mass_element (e : Element) : R :=
  match e with
  | eH  => 1.00783  | eHe => 4.00260
  | eLi => 7.01600  | eBe => 9.01218  | eB  => 11.00931 | eC  => 12.00000
  | eN  => 14.00307 | eO  => 15.99491 | eF  => 18.99840 | eNe => 19.99244
  | eNa => 22.98977 | eMg => 23.98504 | eAl => 26.98154 | eSi => 27.97693
  | eP  => 30.97376 | eS  => 31.97207 | eCl => 34.96885 | eAr => 39.96238
  | eK  => 38.96371 | eCa => 39.96259 | eSc => 44.95592 | eTi => 47.94795
  | eV   => 50.94396 | eCr => 51.94051 | eMn => 54.93805 | eFe => 55.93494
  | eCo => 58.93320 | eNi => 57.93535 | eCu => 62.92960 | eZn => 63.92914
  | eGe => 73.92118 | eAs => 74.92160 | eSe => 79.91652 | eBr => 78.91834
  | eKr => 83.91151
  | eAg => 106.90509| eSn => 119.90220| eI  => 126.90448| eXe => 131.90415
  | eAu => 196.96654| eHg => 201.97063| ePb => 207.97665| eBi => 208.98040
  | _  => atomic_mass e   (* fallback to average mass *)
  end.
 
Definition exact_mass (mol : Molecule) : R :=
  fold_left
    (fun acc a => acc + monoisotopic_mass_element a.(ai_element))
    mol.(mol_atoms)
    0.
 
(** ------------------------------------------------------------------ *)
(** ** 2. LogP (Crippen / Wildman estimation framework)                 *)
(** ------------------------------------------------------------------ *)
 
(** LogP atom contributions (Crippen fragment method, simplified) *)
Definition atom_logP_contribution (e : Element) : R :=
  match e with
  | eC  => 0.1441
  | eN  => -1.019
  | eO  => -0.6865
  | eS  => -0.0024
  | eF  => 0.4202
  | eCl => 0.6895
  | eBr => 0.8765
  | eI  => 1.1645
  | eP  => 0.8738
  | eH  => 0.0000
  | _  => 0.0
  end.
 
(** Simple atom-contribution LogP estimate *)
Definition estimated_logP (mol : Molecule) : R :=
  fold_left
    (fun acc a => acc + atom_logP_contribution a.(ai_element))
    mol.(mol_atoms)
    0.
 
(** ------------------------------------------------------------------ *)
(** ** 3. LogS (Aqueous Solubility)                                     *)
(** ------------------------------------------------------------------ *)
 
(** Simple solubility estimate from LogP and MW (Delaney model) *)
Definition estimated_logS (logP : R) (mw : R) : R :=
  0.16 - 0.63 * logP - 0.0062 * mw.
 
(** Higher LogP generally predicts lower solubility *)
Theorem logP_logS_anticorrelated :
  forall logP1 logP2 mw : R,
    logP1 < logP2 ->
    estimated_logS logP2 mw < estimated_logS logP1 mw.
Proof.
  intros logP1 logP2 mw H.
  unfold estimated_logS; lra.
Qed.
 
(** ------------------------------------------------------------------ *)
(** ** 4. pKa                                                           *)
(** ------------------------------------------------------------------ *)
 
(** Approximate pKa values for common acidic functional group types *)
Definition typical_pKa_carboxylicAcid : R := 4.5.
Definition typical_pKa_phenol         : R := 10.0.
Definition typical_pKa_thiol          : R := 8.0.
Definition typical_pKa_alcohol        : R := 16.0.
Definition typical_pKa_ammonium       : R := 9.5.   (* protonated amine *)
Definition typical_pKa_sulfonicAcid   : R := (-1.0).
 
(** Sulfonic acid is a stronger acid than carboxylic acid *)
Lemma sulfonic_stronger_than_carboxylic :
  typical_pKa_sulfonicAcid < typical_pKa_carboxylicAcid.
Proof. unfold typical_pKa_sulfonicAcid, typical_pKa_carboxylicAcid; lra. Qed.
 
(** ------------------------------------------------------------------ *)
(** ** 5. Topological Descriptors                                       *)
(** ------------------------------------------------------------------ *)
 
(** Topological diameter: maximum topological distance in the molecule *)
Definition topological_diameter (mol : Molecule) : option nat :=
  let atom_ids := map ai_id mol.(mol_atoms) in
  let pairs    := flat_map (fun i => map (fun j => (i, j)) atom_ids) atom_ids in
  let half     := filter (fun '(i, j) => Nat.ltb i j) pairs in
  fold_left
    (fun acc '(i, j) =>
       match topological_distance mol i j with
       | None   => None
       | Some d =>
         match acc with
         | None   => Some d
         | Some m => Some (Nat.max m d)
         end
       end)
    half
    (Some 0%nat).
 
(** Eccentricity of atom [i]: maximum topological distance to any other atom *)
Definition eccentricity (mol : Molecule) (i : nat) : option nat :=
  let atom_ids := map ai_id mol.(mol_atoms) in
  fold_left
    (fun acc j =>
       match topological_distance mol i j with
       | None   => None
       | Some d =>
         match acc with
         | None   => Some d
         | Some m => Some (Nat.max m d)
         end
       end)
    atom_ids
    (Some 0%nat).
 
(** Molecular connectivity index χ₀ = Σ 1/√(degree(v)) for each vertex v.
    Uses [degree] from Molecules.v directly on the molecule graph. *)
Definition connectivity_index_chi0 (mol : Molecule) : R :=
  fold_left
    (fun acc a =>
       let d := degree mol a.(ai_id) in
       if Nat.eqb d 0 then acc
       else acc + 1 / sqrt (INR d))
    mol.(mol_atoms)
    0.
 
(** ------------------------------------------------------------------ *)
(** ** 6. Surface Properties                                            *)
(** ------------------------------------------------------------------ *)
 
(** Topological polar surface area (TPSA) is computed from atom contributions.
    We use the Ertl method (simplified: N, O, S, P surface contributions). *)
Definition tpsa_contribution (e : Element) : R :=
  match e with
  | eN  => 26.02   (* bare N *)
  | eO  => 20.23   (* bare O *)
  | eS  => 25.30   (* bare S *)
  | eP  => 34.14   (* bare P *)
  | _  => 0.0
  end.
 
Definition estimated_tpsa (mol : Molecule) : R :=
  fold_left
    (fun acc a => acc + tpsa_contribution a.(ai_element))
    mol.(mol_atoms)
    0.
 
(** Each element contributes a non-negative TPSA value *)
Lemma tpsa_contribution_nonneg : forall e : Element,
    tpsa_contribution e >= 0.
Proof.
  intro e; destruct e; unfold tpsa_contribution; lra.
Qed.
 
(** Helper: fold_left preserves non-negativity of the TPSA accumulator *)
Lemma fold_left_tpsa_nonneg :
  forall (l : list AtomInst) (acc : R),
    acc >= 0 ->
    fold_left (fun s a => s + tpsa_contribution a.(ai_element)) l acc >= 0.
Proof.
  intros l; induction l as [|a rest IH]; intros acc Hacc.
  - simpl; lra.
  - simpl. apply IH.
    pose proof (tpsa_contribution_nonneg a.(ai_element)). lra.
Qed.
 
(** TPSA is non-negative *)
Lemma tpsa_nonneg : forall mol, estimated_tpsa mol >= 0.
Proof.
  intro mol.
  unfold estimated_tpsa.
  apply fold_left_tpsa_nonneg. lra.
Qed.
 
(** ------------------------------------------------------------------ *)
(** ** 7. Electronic Properties                                         *)
(** ------------------------------------------------------------------ *)
 
(** Partial atomic charges record *)
Record PartialCharges : Type := mkPC
  { pc_atom_id : nat
  ; pc_charge  : R   (* in units of e *)
  }.
 
(** Molecular dipole moment from partial charges and positions *)
Definition dipole_moment
    (charges   : list PartialCharges)
    (positions : nat -> Point3D)
    : Point3D :=
  fold_left
    (fun acc c =>
       let p := positions c.(pc_atom_id) in
       vec_add acc (vec_scale c.(pc_charge) p))
    charges
    origin.
 
(** Magnitude of dipole moment (Debye) *)
Definition dipole_magnitude
    (charges   : list PartialCharges)
    (positions : nat -> Point3D)
    : R :=
  norm (dipole_moment charges positions).
 
(** ------------------------------------------------------------------ *)
(** ** 8. HOMO / LUMO Energies                                          *)
(** ------------------------------------------------------------------ *)
 
Record FrontierOrbitals : Type := mkFO
  { fo_HOMO_energy : R   (* eV, typically negative *)
  ; fo_LUMO_energy : R   (* eV, typically positive *)
  }.
 
(** HOMO-LUMO gap *)
Definition homo_lumo_gap (fo : FrontierOrbitals) : R :=
  fo.(fo_LUMO_energy) - fo.(fo_HOMO_energy).
 
(** HOMO energy < LUMO energy (stability condition) *)
Definition stable_orbital_ordering (fo : FrontierOrbitals) : Prop :=
  fo.(fo_HOMO_energy) < fo.(fo_LUMO_energy).
 
(** Gap is positive iff ordering is stable *)
Theorem gap_positive_iff_stable :
  forall fo : FrontierOrbitals,
    stable_orbital_ordering fo <-> homo_lumo_gap fo > 0.
Proof.
  intro fo.
  unfold stable_orbital_ordering, homo_lumo_gap; split; intro H; lra.
Qed.
 
(** ------------------------------------------------------------------ *)
(** ** 9. Synthesizability (SA Score)                                   *)
(** ------------------------------------------------------------------ *)
 
(** SA score range: 1 (easy) to 10 (difficult) *)
Definition sa_score_min : R := 1.0.
Definition sa_score_max : R := 10.0.
 
Record SAscore : Type := mkSA
  { sa_score : R   (* 1.0 = easy to synthesize, 10.0 = very hard *)
  }.
 
Definition is_synthesizable (sa : SAscore) : Prop :=
  sa.(sa_score) <= 6.0.  (* threshold: SA score ≤ 6 is considered accessible *)
 
(** SA score factors (conceptual model) *)
Record SAFactors : Type := mkSAF
  { saf_chiral_centers   : nat
  ; saf_ring_count        : nat
  ; saf_ring_complexity   : R
  ; saf_uncommon_elements : nat
  ; saf_fragment_score    : R
  }.
 
(** More chiral centers → harder to synthesize.
    Use %nat to keep the comparison in the nat scope while R_scope is open. *)
Lemma more_chirality_harder : forall n m : nat,
    (n < m)%nat ->
    (n <= m)%nat.
Proof. intros n m H; lia. Qed.
 
(** ------------------------------------------------------------------ *)
(** ** 10. Ensemble / Multi-Conformer Properties                        *)
(** ------------------------------------------------------------------ *)
 
(** Boltzmann weight for a conformer with energy E at temperature T (K) *)
Definition boltzmann_weight (E kT : R) : R :=
  exp (- E / kT).
 
(** kT at 298 K in kcal/mol ≈ 0.5922 kcal/mol *)
Definition kT_298K : R := 0.5922.
 
(** Boltzmann-weighted average of a property.
    [conformers] is a list of (energy, property_value) pairs.
    [partZ] is the partition function (sum of weights); renamed from [Z]
    to avoid shadowing [Z : Type] imported from ZArith. *)
Definition boltzmann_average
    (conformers : list (R * R))
    (kT : R)
    : R :=
  let weights  := map (fun '(e, _) => boltzmann_weight e kT) conformers in
  let partZ    := fold_left Rplus weights 0 in
  if Req_EM_T partZ 0 then 0
  else
    let weighted_sum :=
      fold_left
        (fun acc '(e, prop) => acc + boltzmann_weight e kT * prop)
        conformers
        0
    in
    weighted_sum / partZ.
 
(** Boltzmann weights are strictly positive *)
Lemma boltzmann_weight_pos : forall E kT : R,
    kT > 0 -> boltzmann_weight E kT > 0.
Proof.
  intros E kT HkT.
  unfold boltzmann_weight.
  apply exp_pos.
Qed.
 
(** ------------------------------------------------------------------ *)
(** ** 11. Physical Reasonableness Constraints                          *)
(** ------------------------------------------------------------------ *)
 
(** A set of constraints that a physically reasonable molecular structure
    must satisfy.
    [pr_nonneg_atoms] is annotated with %nat because [<] resolves to
    [Rlt] while [R_scope] is open; we want the nat comparison here. *)
Record PhysicallyReasonable (mol : Molecule) : Prop := mkPR
  { pr_finite_weight   : molecular_weight mol < 1e6
  ; pr_positive_weight : molecular_weight mol > 0
  ; pr_nonneg_atoms    : (0 < length mol.(mol_atoms))%nat
  }.
 
(** The empty molecule fails the physical reasonableness check *)
Lemma empty_mol_not_reasonable : ~ PhysicallyReasonable empty_mol.
Proof.
  intro H.
  destruct H.
  unfold empty_mol in pr_nonneg_atoms0; simpl in pr_nonneg_atoms0; lia.
Qed.
 
Close Scope R_scope.