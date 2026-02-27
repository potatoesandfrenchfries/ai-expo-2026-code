(** * DrugLikeness: Lipinski's Rule of Five, Veber's Rules, and Toxicity Filters

    Models:
    - Lipinski's Rule of Five (oral bioavailability)
    - Veber's Rules (oral bioavailability)
    - Ghose filter
    - Egan filter
    - Lead-likeness / fragment-likeness
    - PAINS structural alerts
    - Toxicophores
    - ADMET property ranges
*)

Require Import Stdlib.Reals.Reals.
Require Import Stdlib.Arith.Arith.
Require Import Stdlib.Lists.List.
Require Import Stdlib.micromega.Lra.
Require Import Stdlib.micromega.Lia.
Import ListNotations.
Require Import Chemistry.Atoms.
Require Import Chemistry.Molecules.

Open Scope R_scope.

(** ------------------------------------------------------------------ *)
(** ** 1. Molecular Descriptors Record                                  *)
(** ------------------------------------------------------------------ *)

(** Precomputed molecular descriptors used in drug-likeness filters *)
Record MolDescriptors : Type := mkMolDesc
  { md_mol_weight      : R    (** molecular weight, Da *)
  ; md_logP            : R    (** lipophilicity, log partition coeff *)
  ; md_hbd             : nat  (** H-bond donors (N-H + O-H count) *)
  ; md_hba             : nat  (** H-bond acceptors (N + O count) *)
  ; md_rot_bonds       : nat  (** number of rotatable bonds *)
  ; md_psa             : R    (** polar surface area, Å² *)
  ; md_molar_refract   : R    (** molar refractivity *)
  ; md_logS            : R    (** aqueous solubility, log(mol/L) *)
  ; md_rings           : nat  (** number of rings *)
  ; md_arom_rings      : nat  (** number of aromatic rings *)
  ; md_chiral_centers  : nat  (** number of chiral centers *)
  ; md_atom_count      : nat  (** total atom count (heavy atoms) *)
  }.

(** ------------------------------------------------------------------ *)
(** ** 2. Lipinski's Rule of Five                                       *)
(** ------------------------------------------------------------------ *)

(** Criteria for oral bioavailability (Lipinski 2001) *)
Definition lipinski_mw   (d : MolDescriptors) : Prop := d.(md_mol_weight) <= 500.
Definition lipinski_logP (d : MolDescriptors) : Prop := d.(md_logP) <= 5.
Definition lipinski_hbd  (d : MolDescriptors) : Prop := (d.(md_hbd) <= 5)%nat.
Definition lipinski_hba  (d : MolDescriptors) : Prop := (d.(md_hba) <= 10)%nat.

Definition lipinski_ro5 (d : MolDescriptors) : Prop :=
  lipinski_mw   d /\
  lipinski_logP d /\
  lipinski_hbd  d /\
  lipinski_hba  d.

(** At most one violation is tolerated in the "soft" Rule of Five *)
Definition lipinski_ro5_violations (d : MolDescriptors) : nat :=
  (if Rlt_dec 500 d.(md_mol_weight) then 1 else 0) +
  (if Rlt_dec 5   d.(md_logP)       then 1 else 0) +
  (if Nat.ltb 5   d.(md_hbd)        then 1 else 0) +
  (if Nat.ltb 10  d.(md_hba)        then 1 else 0).

Definition lipinski_soft (d : MolDescriptors) : Prop :=
  (lipinski_ro5_violations d <= 1)%nat.

(** Lipinski rule is strictly more restrictive than soft version *)
Theorem ro5_implies_soft : forall d : MolDescriptors,
    lipinski_ro5 d -> lipinski_soft d.
Proof.
  intros d [Hmw [HlogP [Hhbd Hhba]]].
  unfold lipinski_soft, lipinski_ro5_violations.
  destruct (Rlt_dec 500 (md_mol_weight d)) as [Hmw_gt|Hmw_nlt].
  - exfalso. apply (Rle_not_lt _ _ Hmw). exact Hmw_gt.
  - destruct (Rlt_dec 5 (md_logP d)) as [HlogP_gt|HlogP_nlt].
    + exfalso. apply (Rle_not_lt _ _ HlogP). exact HlogP_gt.
    + assert (Hb : Nat.ltb 5 (md_hbd d) = false).
      { apply Nat.ltb_ge. exact Hhbd. }
      assert (Ha : Nat.ltb 10 (md_hba d) = false).
      { apply Nat.ltb_ge. exact Hhba. }
      rewrite Hb, Ha. simpl. lia.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 3. Veber's Rules (oral bioavailability)                          *)
(** ------------------------------------------------------------------ *)

Definition veber_rot_bonds (d : MolDescriptors) : Prop :=
  (d.(md_rot_bonds) <= 10)%nat.

Definition veber_psa (d : MolDescriptors) : Prop :=
  d.(md_psa) <= 140.

Definition veber_rules (d : MolDescriptors) : Prop :=
  veber_rot_bonds d /\ veber_psa d.

(** ------------------------------------------------------------------ *)
(** ** 4. Ghose Filter                                                  *)
(** ------------------------------------------------------------------ *)

Definition ghose_mw (d : MolDescriptors) : Prop :=
  160 <= d.(md_mol_weight) /\ d.(md_mol_weight) <= 480.

Definition ghose_logP (d : MolDescriptors) : Prop :=
  -0.4 <= d.(md_logP) /\ d.(md_logP) <= 5.6.

Definition ghose_atoms (d : MolDescriptors) : Prop :=
  (20 <= d.(md_atom_count))%nat /\ (d.(md_atom_count) <= 70)%nat.

Definition ghose_molar_refract (d : MolDescriptors) : Prop :=
  40 <= d.(md_molar_refract) /\ d.(md_molar_refract) <= 130.

Definition ghose_filter (d : MolDescriptors) : Prop :=
  ghose_mw d /\ ghose_logP d /\ ghose_atoms d /\ ghose_molar_refract d.

(** ------------------------------------------------------------------ *)
(** ** 5. Egan Filter (absorption)                                      *)
(** ------------------------------------------------------------------ *)

Definition egan_psa (d : MolDescriptors) : Prop :=
  d.(md_psa) <= 131.6.

Definition egan_logP (d : MolDescriptors) : Prop :=
  d.(md_logP) <= 5.88.

Definition egan_filter (d : MolDescriptors) : Prop :=
  egan_psa d /\ egan_logP d.

(** ------------------------------------------------------------------ *)
(** ** 6. Lead-Likeness                                                 *)
(** ------------------------------------------------------------------ *)

(** Oprea lead-like criteria *)
Definition lead_like_mw (d : MolDescriptors) : Prop :=
  d.(md_mol_weight) <= 350.

Definition lead_like_logP (d : MolDescriptors) : Prop :=
  d.(md_logP) <= 3.

Definition lead_like_rings (d : MolDescriptors) : Prop :=
  (d.(md_rings) <= 3)%nat.

Definition lead_like_rot_bonds (d : MolDescriptors) : Prop :=
  (d.(md_rot_bonds) <= 7)%nat.

Definition lead_likeness (d : MolDescriptors) : Prop :=
  lead_like_mw d /\ lead_like_logP d /\
  lead_like_rings d /\ lead_like_rot_bonds d.

(** Lead-like implies lipinski (MW ≤ 350 → MW ≤ 500, logP ≤ 3 → logP ≤ 5) *)
Theorem lead_like_implies_ro5_mw_logP :
  forall d : MolDescriptors,
    (lead_like_mw d -> lipinski_mw d) /\
    (lead_like_logP d -> lipinski_logP d).
Proof.
  intros d; split; intro H.
  - unfold lead_like_mw, lipinski_mw in *; lra.
  - unfold lead_like_logP, lipinski_logP in *; lra.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 7. Fragment-Likeness (Rule of Three)                             *)
(** ------------------------------------------------------------------ *)

Definition fragment_like_mw (d : MolDescriptors) : Prop :=
  d.(md_mol_weight) <= 300.

Definition fragment_like_logP (d : MolDescriptors) : Prop :=
  d.(md_logP) <= 3.

Definition fragment_like_hbd (d : MolDescriptors) : Prop :=
  (d.(md_hbd) <= 3)%nat.

Definition fragment_like_hba (d : MolDescriptors) : Prop :=
  (d.(md_hba) <= 3)%nat.

Definition fragment_like_rot_bonds (d : MolDescriptors) : Prop :=
  (d.(md_rot_bonds) <= 3)%nat.

Definition fragment_rule_of_three (d : MolDescriptors) : Prop :=
  fragment_like_mw d /\ fragment_like_logP d /\
  fragment_like_hbd d /\ fragment_like_hba d /\
  fragment_like_rot_bonds d.

(** Fragment-like plus a ring bound implies lead-like *)
Theorem fragment_implies_lead_like :
  forall d : MolDescriptors,
    fragment_rule_of_three d ->
    lead_like_rings d ->
    lead_likeness d.
Proof.
  intros d [Hmw [HlogP [_ [_ Hrot]]]] Hrings.
  unfold lead_likeness, lead_like_mw, lead_like_logP,
         lead_like_rings, lead_like_rot_bonds in *.
  repeat split.
  - eapply Rle_trans; [exact Hmw | lra].
  - eapply Rle_trans; [exact HlogP | lra].
  - exact Hrings.
  - eapply Nat.le_trans; [exact Hrot | lia].
Qed.

(** ------------------------------------------------------------------ *)
(** ** 8. PAINS Structural Alerts                                       *)
(** ------------------------------------------------------------------ *)

(** Pan-Assay Interference Compounds: structural features that cause
    false positives in HTS assays.  We model each alert as a boolean
    flag derived from substructure analysis (the actual substructure
    matching is left abstract). *)

Inductive PAINSAlert : Type :=
  | ReactiveElectrophile
  | MichaelAcceptor
  | AldehydeGroup
  | AcylHalide
  | EpoxideGroup
  | PeroxideGroup
  | QuinoneGroup
  | FreeThiol
  | NitroGroup
  | HydrazineGroup
  | IsocyanateGroup
  | RhodanineCore
  | SalicylamideCore.

(** A molecule passes the PAINS filter if it has no alerts *)
Definition passes_pains (alerts : list PAINSAlert) : Prop :=
  alerts = [].

(** ------------------------------------------------------------------ *)
(** ** 9. Toxicophores                                                  *)
(** ------------------------------------------------------------------ *)

Inductive Toxicophore : Type :=
  | AromaticAmine
  | NitroAromatic
  | PolyhalogenatedCompound
  | HeavyMetal
  | AzoCompound
  | ThioureaGroup
  | MichaelAcceptorTox
  | Carcinogen_PAH   (** polycyclic aromatic hydrocarbons *)
  | Mutagen_Ames.    (** Ames test positive structural features *)

Definition is_toxic (tox : list Toxicophore) : Prop :=
  exists t, In t tox.

(** A clean molecule has no toxicophores *)
Definition no_toxicophores (tox : list Toxicophore) : Prop :=
  tox = [].

(** ------------------------------------------------------------------ *)
(** ** 10. ADMET Property Ranges                                        *)
(** ------------------------------------------------------------------ *)

Record ADMETProfile : Type := mkADMET
  { adm_solubility_logS    : R    (** aqueous solubility, log(mol/L); > -4 good *)
  ; adm_caco2_perm         : R    (** Caco-2 permeability, nm/s; > 20 good *)
  ; adm_bbb                : bool (** blood-brain barrier penetration *)
  ; adm_pgp_substrate      : bool (** P-glycoprotein substrate (efflux) *)
  ; adm_cyp3a4_inhibitor   : bool (** CYP3A4 inhibition *)
  ; adm_herg_inhibitor      : bool (** hERG channel inhibition (QT prolongation) *)
  ; adm_hepatotoxic         : bool (** hepatotoxicity *)
  ; adm_mutagenic_ames      : bool (** Ames mutagenicity *)
  ; adm_carcinogenic        : bool
  }.

(** Good ADMET profile: no major liabilities *)
Definition good_admet (adm : ADMETProfile) : Prop :=
  adm.(adm_solubility_logS) >= -4  /\
  adm.(adm_caco2_perm)     >= 20   /\
  adm.(adm_herg_inhibitor) = false /\
  adm.(adm_hepatotoxic)    = false /\
  adm.(adm_mutagenic_ames) = false /\
  adm.(adm_carcinogenic)   = false.

(** ------------------------------------------------------------------ *)
(** ** 11. Combined Drug Candidate Filter                               *)
(** ------------------------------------------------------------------ *)

Definition is_drug_candidate
    (d    : MolDescriptors)
    (alerts: list PAINSAlert)
    (tox  : list Toxicophore)
    (adm  : ADMETProfile)
    : Prop :=
  lipinski_ro5 d     /\
  veber_rules  d     /\
  passes_pains alerts/\
  no_toxicophores tox/\
  good_admet adm.

Close Scope R_scope.
