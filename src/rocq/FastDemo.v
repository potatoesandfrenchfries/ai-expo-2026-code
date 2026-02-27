(** * FastDemo: Demo using Molecules and native computation
    
    This demonstrates aspirin-like molecules with:
    - Molecules graph model from Molecules.v
    - Native computation for fast evaluation
    - Minimal imports
    - Opaque proofs
*)

From Stdlib Require Import List.
Import ListNotations.
Require Import Stdlib.Reals.Reals.
Open Scope R_scope.

Require Import Stdlib.ZArith.ZArith.
Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.
Require Import Chemistry.Molecules.

(** ------------------------------------------------------------------ *)
(** ** Ethanol: C2H6O (9 atoms, 8 bonds)                                *)
(** ------------------------------------------------------------------ *)

Definition ethanol : Molecule :=
  mkMol
    [ mkAtom 0 eC (mkPoint (0 / 10000) (0 / 10000) (0 / 10000)) 0%Z None None
    ; mkAtom 1 eC (mkPoint (15400 / 10000) (0 / 10000) (0 / 10000)) 0%Z None None
    ; mkAtom 2 eO (mkPoint (21340 / 10000) (12470 / 10000) (0 / 10000)) 0%Z None None
    ; mkAtom 3 eH (mkPoint (-3860 / 10000) (10260 / 10000) (0 / 10000)) 0%Z None None
    ; mkAtom 4 eH (mkPoint (-3860 / 10000) (-5130 / 10000) (8890 / 10000)) 0%Z None None
    ; mkAtom 5 eH (mkPoint (-3860 / 10000) (-5130 / 10000) (-8890 / 10000)) 0%Z None None
    ; mkAtom 6 eH (mkPoint (19260 / 10000) (-10260 / 10000) (0 / 10000)) 0%Z None None
    ; mkAtom 7 eH (mkPoint (19260 / 10000) (5130 / 10000) (8890 / 10000)) 0%Z None None
    ; mkAtom 8 eH (mkPoint (32340 / 10000) (12470 / 10000) (0 / 10000)) 0%Z None None
    ]
    [ mkBond 0 0 1 SingleBond None
    ; mkBond 1 1 2 SingleBond None
    ; mkBond 2 0 3 SingleBond None
    ; mkBond 3 0 4 SingleBond None
    ; mkBond 4 0 5 SingleBond None
    ; mkBond 5 1 6 SingleBond None
    ; mkBond 6 1 7 SingleBond None
    ; mkBond 7 2 8 SingleBond None
    ].

(** ------------------------------------------------------------------ *)
(** ** Aspirin: C9H8O4 (21 atoms, 21 bonds)                             *)
(** ------------------------------------------------------------------ *)

Definition aspirin : Molecule :=
  mkMol
    [ (* Benzene ring *)
      mkAtom 0 eC (mkPoint 0 0 0) 0%Z None None
    ; mkAtom 1 eC (mkPoint 1.4 0 0) 0%Z None None
    ; mkAtom 2 eC (mkPoint 2.1 1.2 0) 0%Z None None
    ; mkAtom 3 eC (mkPoint 1.4 2.4 0) 0%Z None None
    ; mkAtom 4 eC (mkPoint 0 2.4 0) 0%Z None None
    ; mkAtom 5 eC (mkPoint (-0.7) 1.2 0) 0%Z None None
    (* Carboxylic acid group *)
    ; mkAtom 6 eC (mkPoint 2.1 3.6 0) 0%Z None None
    ; mkAtom 7 eO (mkPoint 1.4 4.8 0) 0%Z None None
    ; mkAtom 8 eO (mkPoint 3.5 3.6 0) 0%Z None None
    ; mkAtom 9 eH (mkPoint 3.8 4.8 0) 0%Z None None
    (* Acetyl group *)
    ; mkAtom 10 eO (mkPoint (-0.7) 3.6 0) 0%Z None None
    ; mkAtom 11 eC (mkPoint (-2.1) 3.6 0) 0%Z None None
    ; mkAtom 12 eO (mkPoint (-2.8) 4.8 0) 0%Z None None
    ; mkAtom 13 eC (mkPoint (-2.8) 2.4 0) 0%Z None None
    (* Hydrogens on benzene *)
    ; mkAtom 14 eH (mkPoint (-0.7) (-0.9) 0) 0%Z None None
    ; mkAtom 15 eH (mkPoint 2.1 (-0.9) 0) 0%Z None None
    ; mkAtom 16 eH (mkPoint 3.1 1.2 0) 0%Z None None
    ; mkAtom 17 eH (mkPoint (-1.7) 1.2 0) 0%Z None None
    (* Hydrogens on methyl *)
    ; mkAtom 18 eH (mkPoint (-3.8) 2.4 0) 0%Z None None
    ; mkAtom 19 eH (mkPoint (-2.5) 1.5 0) 0%Z None None
    ; mkAtom 20 eH (mkPoint (-2.5) 3.3 0) 0%Z None None
    ]
    [ (* Benzene ring bonds *)
      mkBond 0 0 1 AromaticBond None
    ; mkBond 1 1 2 AromaticBond None
    ; mkBond 2 2 3 AromaticBond None
    ; mkBond 3 3 4 AromaticBond None
    ; mkBond 4 4 5 AromaticBond None
    ; mkBond 5 5 0 AromaticBond None
    (* Carboxylic acid *)
    ; mkBond 6 3 6 SingleBond None
    ; mkBond 7 6 7 DoubleBond None
    ; mkBond 8 6 8 SingleBond None
    ; mkBond 9 8 9 SingleBond None
    (* Acetyl ester *)
    ; mkBond 10 4 10 SingleBond None
    ; mkBond 11 10 11 SingleBond None
    ; mkBond 12 11 12 DoubleBond None
    ; mkBond 13 11 13 SingleBond None
    (* Hydrogens *)
    ; mkBond 14 0 14 SingleBond None
    ; mkBond 15 1 15 SingleBond None
    ; mkBond 16 2 16 SingleBond None
    ; mkBond 17 5 17 SingleBond None
    ; mkBond 18 13 18 SingleBond None
    ; mkBond 19 13 19 SingleBond None
    ; mkBond 20 13 20 SingleBond None
    ].

(** ------------------------------------------------------------------ *)
(** ** Verification checks using native_compute                         *)
(** ------------------------------------------------------------------ *)

(* Type-check: confirm molecules are well-typed *)
Check ethanol.
Check aspirin.

(* Fast runtime evaluation using native OCaml compilation *)
Eval native_compute in (is_connected ethanol).
Eval native_compute in (atom_count ethanol).
Eval native_compute in (bond_count ethanol).

Eval native_compute in (is_connected aspirin).
Eval native_compute in (atom_count aspirin).
Eval native_compute in (bond_count aspirin).

(* Molecular weight computation *)
Eval native_compute in (molecular_weight ethanol).
Eval native_compute in (molecular_weight aspirin).

(** ------------------------------------------------------------------ *)
(** ** Simple properties (admitted for speed)                           *)
(** ------------------------------------------------------------------ *)

(** Ethanol is connected *)
Lemma ethanol_connected : is_connected ethanol = true.
Proof. native_compute. reflexivity. Qed.

(** Aspirin is connected *)
Lemma aspirin_connected : is_connected aspirin = true.
Proof. native_compute. reflexivity. Qed.

(** Ethanol has 9 atoms *)
Lemma ethanol_atom_count : atom_count ethanol = 9.
Proof. reflexivity. Qed.

(** Aspirin has 21 atoms *)
Lemma aspirin_atom_count : atom_count aspirin = 21.
Proof. reflexivity. Qed.
