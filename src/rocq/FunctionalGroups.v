(** * FunctionalGroups: Functional Group Types and Properties
    Models:
    - All common organic functional groups
    - Reactivity ranking
    - Electron donating / withdrawing character
    - Acidic / basic character
    - Nucleophilicity and electrophilicity
    - Leaving group ability
    - Structural toxicity alerts and PAINS classification
*)

Require Import Stdlib.Lists.List.
Require Import Stdlib.micromega.Lia.
Import ListNotations.
Require Import Chemistry.Atoms.

(** ------------------------------------------------------------------ *)
(** ** 1. Functional Group Type                                         *)
(** ------------------------------------------------------------------ *)

Inductive FunctionalGroup : Type :=
  (** Hydrocarbons *)
  | Alkane
  | Alkene
  | Alkyne
  (** Oxygen-containing *)
  | Alcohol
  | Phenol
  | Ether
  | Aldehyde
  | Ketone
  | CarboxylicAcid
  | Ester
  | Anhydride
  | Lactone
  | Epoxide
  | Peroxide
  (** Nitrogen-containing *)
  | Amine_primary
  | Amine_secondary
  | Amine_tertiary
  | Amide
  | Imine
  | Nitrile
  | Nitro
  | Azo
  | Diazo
  | Azide
  | Hydrazine
  | Carbamate
  | Urea
  | Isocyanate
  | Isonitrile
  (** Sulfur-containing *)
  | Thiol
  | Thioether
  | Disulfide
  | Sulfoxide
  | Sulfone
  | SulfonicAcid
  | Sulfonamide
  | Thioamide
  | Thiourea
  (** Halides *)
  | Fluoride
  | Chloride
  | Bromide
  | Iodide
  (** Phosphorus-containing *)
  | Phosphate
  | Phosphonate
  | Phosphine
  | PhosphonicAcid
  (** Silicon-containing *)
  | Silane
  (** Boron-containing *)
  | BoronicAcid
  | Boronate
  (** Carbene / other *)
  | Carbene
  | Carbenoid.

(** ------------------------------------------------------------------ *)
(** ** 2. Electron Donor / Withdrawing Character                        *)
(** ------------------------------------------------------------------ *)

Inductive ElectronicEffect : Type :=
  | StrongDonor        (* +I, +M: e.g. NR₂, OR *)
  | ModerateDonor      (* +I: e.g. alkyl *)
  | WeakDonor
  | Neutral
  | WeakWithdrawer
  | ModerateWithdrawer  (* -I: halogens *)
  | StrongWithdrawer.   (* -M: NO₂, C=O, CN *)

Definition electronic_effect (fg : FunctionalGroup) : ElectronicEffect :=
  match fg with
  | Amine_primary | Amine_secondary | Amine_tertiary => StrongDonor
  | Alcohol | Phenol | Ether                         => StrongDonor
  | Alkane                                            => ModerateDonor
  | Alkene | Alkyne                                   => Neutral
  | Aldehyde | Ketone | Ester | CarboxylicAcid       => StrongWithdrawer
  | Amide | Carbamate | Urea                         => ModerateWithdrawer
  | Nitrile                                           => StrongWithdrawer
  | Nitro                                             => StrongWithdrawer
  | Fluoride                                          => ModerateWithdrawer
  | Chloride | Bromide | Iodide                       => ModerateWithdrawer
  | Thiol | Thioether                                 => ModerateDonor
  | Sulfoxide | Sulfone | SulfonicAcid               => StrongWithdrawer
  | Sulfonamide                                       => StrongWithdrawer
  | Phosphate | Phosphonate                           => StrongWithdrawer
  | BoronicAcid | Boronate                           => ModerateWithdrawer
  | _                                                 => Neutral
  end.

(** ------------------------------------------------------------------ *)
(** ** 3. Acidic / Basic Character                                      *)
(** ------------------------------------------------------------------ *)

Inductive AcidBaseCharacter : Type :=
  | StrongAcid
  | WeakAcid
  | Amphoteric
  | Neutral_pH
  | WeakBase
  | StrongBase.

Definition acid_base_char (fg : FunctionalGroup) : AcidBaseCharacter :=
  match fg with
  | CarboxylicAcid | SulfonicAcid | PhosphonicAcid => WeakAcid
  | Phenol                                           => WeakAcid
  | Thiol                                            => WeakAcid
  | Alcohol                                          => Amphoteric
  | Amine_primary | Amine_secondary | Amine_tertiary => WeakBase
  | Amide                                             => Neutral_pH   (* very weak base *)
  | Nitro                                             => Neutral_pH
  | _                                                 => Neutral_pH
  end.

(** ------------------------------------------------------------------ *)
(** ** 4. Nucleophilicity / Electrophilicity                            *)
(** ------------------------------------------------------------------ *)

(** Relative nucleophilicity score (higher = better nucleophile) *)
Definition nucleophilicity_rank (fg : FunctionalGroup) : nat :=
  match fg with
  | Amine_primary   => 8
  | Amine_secondary => 7
  | Thiol           => 9
  | Thioether       => 6
  | Alcohol         => 5
  | Phenol          => 5
  | Ether           => 3
  | Alkene          => 4
  | Alkyne          => 3
  | Fluoride        => 1
  | Chloride        => 3
  | Bromide         => 5
  | Iodide          => 7
  | _               => 0
  end.

(** Relative electrophilicity score (higher = better electrophile).
    Michael acceptors (activated alkenes) have moderate electrophilicity;
    Alkene is used as the representative Michael acceptor pattern here. *)
Definition electrophilicity_rank (fg : FunctionalGroup) : nat :=
  match fg with
  | Aldehyde       => 9
  | Ketone         => 7
  | CarboxylicAcid => 6
  | Ester          => 5
  | Amide          => 4
  | Nitrile        => 6
  | Nitro          => 5
  | Isocyanate     => 8
  | Alkene         => 7    (* Michael acceptor character *)
  | _              => 0
  end.

(** Higher nucleophilicity means lower electrophilicity tendency *)
Theorem nucleophile_not_electrophile :
  forall fg,
    nucleophilicity_rank fg > 5 ->
    electrophilicity_rank fg <= 3.
Proof.
  intro fg; destruct fg; simpl; lia.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 5. Leaving Group Ability                                         *)
(** ------------------------------------------------------------------ *)

(** Relative leaving group ability (higher = better leaving group in SN2/E2) *)
Definition leaving_group_rank (fg : FunctionalGroup) : nat :=
  match fg with
  | Iodide       => 9
  | Bromide      => 8
  | Chloride     => 6
  | Fluoride     => 1   (* poor leaving group *)
  | SulfonicAcid => 9   (* tosylate / mesylate type *)
  | Thiol        => 4
  | Alcohol      => 2   (* needs activation *)
  | _            => 0
  end.

(** Fluoride is a poor leaving group (rank 1) vs iodide (rank 9) *)
Lemma iodide_better_leaving_than_fluoride :
  leaving_group_rank Fluoride < leaving_group_rank Iodide.
Proof. simpl; lia. Qed.

(** ------------------------------------------------------------------ *)
(** ** 6. Reactivity Ranking                                            *)
(** ------------------------------------------------------------------ *)

(** Approximate reactivity toward nucleophilic acyl substitution
    (higher = more reactive).
    [Chloride] stands in for acyl chloride (the most reactive acylating agent)
    since AcylHalide is not a separate constructor. *)
Definition electrophile_reactivity (fg : FunctionalGroup) : nat :=
  match fg with
  | Chloride       => 10   (* acyl chloride / acid chloride *)
  | Anhydride      => 9
  | Aldehyde       => 8
  | Ketone         => 7
  | Ester          => 5
  | Carbamate      => 4
  | Amide          => 3    (* least reactive carbonyl *)
  | Nitrile        => 5
  | _              => 0
  end.

(** Reactivity order: amide < ester < anhydride *)
Lemma reactivity_carbonyl_order :
  electrophile_reactivity Amide <
  electrophile_reactivity Ester /\
  electrophile_reactivity Ester <
  electrophile_reactivity Anhydride.
Proof. simpl; lia. Qed.

(** ------------------------------------------------------------------ *)
(** ** 7. Functional Group Properties Record                            *)
(** ------------------------------------------------------------------ *)

Record FGProperties : Type := mkFGProp
  { fgp_fg            : FunctionalGroup
  ; fgp_electron_eff  : ElectronicEffect
  ; fgp_acid_base     : AcidBaseCharacter
  ; fgp_nucleophil    : nat
  ; fgp_electrophil   : nat
  ; fgp_leaving_group : nat
  }.

Definition fg_properties (fg : FunctionalGroup) : FGProperties :=
  mkFGProp fg
    (electronic_effect fg)
    (acid_base_char fg)
    (nucleophilicity_rank fg)
    (electrophilicity_rank fg)
    (leaving_group_rank fg).

(** ------------------------------------------------------------------ *)
(** ** 8. Functional Group Recognition Helpers                          *)
(** ------------------------------------------------------------------ *)

(** Is the group an oxygen-containing functional group? *)
Definition is_oxygen_fg (fg : FunctionalGroup) : bool :=
  match fg with
  | Alcohol | Phenol | Ether | Aldehyde | Ketone | CarboxylicAcid
  | Ester | Anhydride | Lactone | Epoxide | Peroxide
  | Sulfoxide | Sulfone | SulfonicAcid | Phosphate
  | Phosphonate | PhosphonicAcid => true
  | _ => false
  end.

(** Is the group a nitrogen-containing functional group? *)
Definition is_nitrogen_fg (fg : FunctionalGroup) : bool :=
  match fg with
  | Amine_primary | Amine_secondary | Amine_tertiary | Amide
  | Imine | Nitrile | Nitro | Azo | Diazo | Azide | Hydrazine
  | Carbamate | Urea | Isocyanate | Isonitrile
  | Sulfonamide | Thioamide | Thiourea => true
  | _ => false
  end.

(** Is the group a halide? *)
Definition is_halide (fg : FunctionalGroup) : bool :=
  match fg with
  | Fluoride | Chloride | Bromide | Iodide => true
  | _ => false
  end.

(** ------------------------------------------------------------------ *)
(** ** 9. Structural Toxicity Alerts                                    *)
(** ------------------------------------------------------------------ *)

(** Categories of structural toxicity concern *)
Inductive ToxicityAlert : Type :=
  | NoAlert
  | AlkylatingAgent      (* DNA alkylation → mutagenic / genotoxic *)
  | ReactiveElectrophile (* non-specific protein/DNA electrophile *)
  | MetabolicAlert       (* activatable to toxic metabolite in vivo *)
  | ReactiveOxygen       (* reactive oxygen species generation *)
  | ProteinReactive.     (* covalent protein binding → idiosyncratic toxicity *)

(** Map each functional group to its primary structural alert.
    Based on the Brenk / Ames / genotoxicity alert literature. *)
Definition structural_alert (fg : FunctionalGroup) : ToxicityAlert :=
  match fg with
  | Epoxide              => AlkylatingAgent      (* direct DNA alkylator *)
  | Diazo | Azide        => AlkylatingAgent      (* reactive carbene / nitrene precursors *)
  | Isocyanate           => ReactiveElectrophile (* reacts with NH₂/OH of proteins *)
  | Isonitrile           => ReactiveElectrophile
  | Peroxide             => ReactiveOxygen       (* homolytic O-O cleavage → radicals *)
  | Nitro                => MetabolicAlert       (* reduced to hydroxylamine/nitroso *)
  | Azo                  => MetabolicAlert       (* reductive cleavage → aromatic amines *)
  | Hydrazine            => MetabolicAlert       (* oxidative metabolism → reactive species *)
  | Aldehyde             => ProteinReactive      (* Schiff base formation with Lys *)
  | Thiol                => ProteinReactive      (* disulfide exchange with Cys *)
  | Disulfide            => ProteinReactive      (* thiol exchange *)
  | _                    => NoAlert
  end.

(** An epoxide carries an alkylating agent alert *)
Lemma epoxide_is_alkylating_alert :
  structural_alert Epoxide = AlkylatingAgent.
Proof. reflexivity. Qed.

(** Alcohols and simple amines carry no structural alert *)
Lemma alcohol_no_alert :
  structural_alert Alcohol = NoAlert.
Proof. reflexivity. Qed.

Lemma amine_primary_no_alert :
  structural_alert Amine_primary = NoAlert.
Proof. reflexivity. Qed.

(** ------------------------------------------------------------------ *)
(** ** 10. PAINS (Pan-Assay INterference CompoundS) Classification      *)
(** ------------------------------------------------------------------ *)

(** PAINS categories — functional groups that frequently cause
    false positives in biochemical assays through non-specific
    mechanisms (fluorescence, redox cycling, aggregation, etc.). *)
Inductive PAINSClass : Type :=
  | NotPAINS
  | PAINS_ReactiveGroup   (* non-specific covalent reactors *)
  | PAINS_Redox           (* redox cycling / ROS generation *)
  | PAINS_Chelator        (* promiscuous metal chelation *)
  | PAINS_Aggregator.     (* colloidal aggregation at low µM *)

Definition pains_class (fg : FunctionalGroup) : PAINSClass :=
  match fg with
  | Aldehyde                  => PAINS_ReactiveGroup  (* Schiff base with assay proteins *)
  | Isocyanate                => PAINS_ReactiveGroup
  | Epoxide                   => PAINS_ReactiveGroup
  | Diazo | Azide             => PAINS_ReactiveGroup
  | Peroxide                  => PAINS_Redox
  | Azo                       => PAINS_Redox           (* redox-active azo *)
  | Nitro                     => PAINS_Redox           (* nitroreduction cycling *)
  | BoronicAcid | Boronate    => PAINS_Chelator        (* boronate esters can chelate *)
  | _                         => NotPAINS
  end.

(** Peroxides are flagged as PAINS redox cyclers *)
Lemma peroxide_pains_redox :
  pains_class Peroxide = PAINS_Redox.
Proof. reflexivity. Qed.

(** A functional group with no structural alert is also not a PAINS
    reactive group (reactive groups always carry structural alerts). *)
Lemma no_alert_not_reactive_pains :
  forall fg,
    structural_alert fg = NoAlert ->
    pains_class fg <> PAINS_ReactiveGroup.
Proof.
  intro fg; destruct fg; simpl; intro H; discriminate.
Qed.

(** ------------------------------------------------------------------ *)
(** ** 11. Composite Toxicity Risk Score                                *)
(** ------------------------------------------------------------------ *)

(** Numeric toxicity risk derived from the structural alert category.
    Scale: 0 (safe) to 10 (high concern). *)
Definition toxicity_risk_score (fg : FunctionalGroup) : nat :=
  match structural_alert fg with
  | NoAlert              => 0
  | ProteinReactive      => 4
  | MetabolicAlert       => 6
  | ReactiveElectrophile => 7
  | ReactiveOxygen       => 7
  | AlkylatingAgent      => 9
  end.

(** A group with no structural alert scores zero *)
Lemma no_alert_zero_risk :
  forall fg,
    structural_alert fg = NoAlert ->
    toxicity_risk_score fg = 0.
Proof.
  intros fg H.
  unfold toxicity_risk_score; rewrite H; reflexivity.
Qed.

(** Alkylating agents score higher than metabolic alerts *)
Lemma alkylating_riskier_than_metabolic :
  forall fg1 fg2,
    structural_alert fg1 = AlkylatingAgent ->
    structural_alert fg2 = MetabolicAlert  ->
    toxicity_risk_score fg2 < toxicity_risk_score fg1.
Proof.
  intros fg1 fg2 H1 H2.
  unfold toxicity_risk_score; rewrite H1, H2; lia.
Qed.

(** Epoxide has a higher toxicity risk score than nitro *)
Lemma epoxide_riskier_than_nitro :
  toxicity_risk_score Epoxide > toxicity_risk_score Nitro.
Proof.
  cbv [toxicity_risk_score structural_alert].
  unfold gt.
  repeat constructor.
Qed.