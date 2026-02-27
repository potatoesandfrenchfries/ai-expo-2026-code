From Stdlib Require Import List.
Import ListNotations.
Require Import Stdlib.Reals.Reals.
Open Scope R_scope.

Require Import Stdlib.ZArith.ZArith.
Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.
Require Import Chemistry.Molecules.

Require Import Chemistry.MolecularProperties.
Require Import Chemistry.HydrogenBonding.
Require Import Chemistry.FunctionalGroups.
Require Import Chemistry.MathProperties.
Require Import Chemistry.Valency.
Require Import Chemistry.Aromaticity.
Require Import Chemistry.Conformational.

Definition demo_molecule : Molecule :=
  mkMol
    [ mkAtom 0 eC (mkPoint (-57088 / 10000) (17830 / 10000) (-6281 / 10000)) 0%Z None None
    ; mkAtom 1 eN (mkPoint (-50537 / 10000) (4348 / 10000) (-7718 / 10000)) 1%Z None None
    ; mkAtom 2 eC (mkPoint (-44740 / 10000) (-386 / 10000) (5473 / 10000)) 0%Z None None
    ; mkAtom 3 eC (mkPoint (-34951 / 10000) (-11694 / 10000) (3019 / 10000)) 0%Z None None
    ; mkAtom 4 eC (mkPoint (-22940 / 10000) (-7011 / 10000) (-5407 / 10000)) 0%Z None None
    ; mkAtom 5 eN (mkPoint (-11424 / 10000) (-3532 / 10000) (2854 / 10000)) 0%Z None None
    ; mkAtom 6 eC (mkPoint (-4213 / 10000) (-13427 / 10000) (9383 / 10000)) 0%Z None None
    ; mkAtom 7 eO (mkPoint (-8371 / 10000) (-24942 / 10000) (10333 / 10000)) 0%Z None None
    ; mkAtom 8 eC (mkPoint (8522 / 10000) (-8686 / 10000) (16087 / 10000)) 0%Z None None
    ; mkAtom 9 eC (mkPoint (17672 / 10000) (-20146 / 10000) (20177 / 10000)) 0%Z None None
    ; mkAtom 10 eC (mkPoint (31414 / 10000) (-13792 / 10000) (20373 / 10000)) 0%Z None None
    ; mkAtom 11 eN (mkPoint (30802 / 10000) (-3583 / 10000) (10140 / 10000)) 0%Z None None
    ; mkAtom 12 eC (mkPoint (42437 / 10000) (2726 / 10000) (4901 / 10000)) 0%Z None None
    ; mkAtom 13 eC (mkPoint (41725 / 10000) (12163 / 10000) (-5425 / 10000)) 0%Z None None
    ; mkAtom 14 eC (mkPoint (53257 / 10000) (18297 / 10000) (-10505 / 10000)) 0%Z None None
    ; mkAtom 15 eC (mkPoint (65799 / 10000) (15107 / 10000) (-5353 / 10000)) 0%Z None None
    ; mkAtom 16 eCl (mkPoint (79573 / 10000) (23029 / 10000) (-12075 / 10000)) 0%Z None None
    ; mkAtom 17 eC (mkPoint (66749 / 10000) (5738 / 10000) (4928 / 10000)) 0%Z None None
    ; mkAtom 18 eCl (mkPoint (81887 / 10000) (1206 / 10000) (11940 / 10000)) 0%Z None None
    ; mkAtom 19 eC (mkPoint (55163 / 10000) (-367 / 10000) (9981 / 10000)) 0%Z None None
    ; mkAtom 20 eC (mkPoint (17485 / 10000) (-294 / 10000) (7266 / 10000)) 0%Z None None
    ; mkAtom 21 eO (mkPoint (13043 / 10000) (8099 / 10000) (-479 / 10000)) 0%Z None None
    ; mkAtom 22 eC (mkPoint (-26586 / 10000) (4805 / 10000) (-14508 / 10000)) 0%Z None None
    ; mkAtom 23 eC (mkPoint (-40903 / 10000) (3722 / 10000) (-19400 / 10000)) 0%Z None None
    ; mkAtom 24 eH (mkPoint (-62401 / 10000) (20010 / 10000) (-15581 / 10000)) 0%Z None None
    ; mkAtom 25 eH (mkPoint (-49327 / 10000) (25279 / 10000) (-4342 / 10000)) 0%Z None None
    ; mkAtom 26 eH (mkPoint (-64181 / 10000) (17278 / 10000) (2017 / 10000)) 0%Z None None
    ; mkAtom 27 eH (mkPoint (-58164 / 10000) (-2189 / 10000) (-10001 / 10000)) 0%Z None None
    ; mkAtom 28 eH (mkPoint (-39863 / 10000) (8152 / 10000) (10317 / 10000)) 0%Z None None
    ; mkAtom 29 eH (mkPoint (-53193 / 10000) (-3657 / 10000) (11617 / 10000)) 0%Z None None
    ; mkAtom 30 eH (mkPoint (-31827 / 10000) (-15868 / 10000) (12646 / 10000)) 0%Z None None
    ; mkAtom 31 eH (mkPoint (-40129 / 10000) (-19771 / 10000) (-2321 / 10000)) 0%Z None None
    ; mkAtom 32 eH (mkPoint (-19787 / 10000) (-15447 / 10000) (-11694 / 10000)) 0%Z None None
    ; mkAtom 33 eH (mkPoint (-5580 / 10000) (4292 / 10000) (-257 / 10000)) 0%Z None None
    ; mkAtom 34 eH (mkPoint (5807 / 10000) (-2664 / 10000) (24844 / 10000)) 0%Z None None
    ; mkAtom 35 eH (mkPoint (17431 / 10000) (-28215 / 10000) (12736 / 10000)) 0%Z None None
    ; mkAtom 36 eH (mkPoint (14960 / 10000) (-24494 / 10000) (29848 / 10000)) 0%Z None None
    ; mkAtom 37 eH (mkPoint (33434 / 10000) (-8982 / 10000) (30011 / 10000)) 0%Z None None
    ; mkAtom 38 eH (mkPoint (39083 / 10000) (-21287 / 10000) (18194 / 10000)) 0%Z None None
    ; mkAtom 39 eH (mkPoint (32278 / 10000) (15009 / 10000) (-9948 / 10000)) 0%Z None None
    ; mkAtom 40 eH (mkPoint (52316 / 10000) (25576 / 10000) (-18534 / 10000)) 0%Z None None
    ; mkAtom 41 eH (mkPoint (56399 / 10000) (-7506 / 10000) (18071 / 10000)) 0%Z None None
    ; mkAtom 42 eH (mkPoint (-25218 / 10000) (14348 / 10000) (-9250 / 10000)) 0%Z None None
    ; mkAtom 43 eH (mkPoint (-19712 / 10000) (5048 / 10000) (-23052 / 10000)) 0%Z None None
    ; mkAtom 44 eH (mkPoint (-42669 / 10000) (-5873 / 10000) (-24392 / 10000)) 0%Z None None
    ; mkAtom 45 eH (mkPoint (-43432 / 10000) (11748 / 10000) (-26404 / 10000)) 0%Z None None
    ]
    [ mkBond 0 0 1 SingleBond None
    ; mkBond 1 1 2 SingleBond None
    ; mkBond 2 2 3 SingleBond None
    ; mkBond 3 3 4 SingleBond None
    ; mkBond 4 4 5 SingleBond None
    ; mkBond 5 5 6 SingleBond None
    ; mkBond 6 6 7 DoubleBond None
    ; mkBond 7 6 8 SingleBond None
    ; mkBond 8 8 9 SingleBond None
    ; mkBond 9 9 10 SingleBond None
    ; mkBond 10 10 11 SingleBond None
    ; mkBond 11 11 12 SingleBond None
    ; mkBond 12 12 13 AromaticBond None
    ; mkBond 13 13 14 AromaticBond None
    ; mkBond 14 14 15 AromaticBond None
    ; mkBond 15 15 16 SingleBond None
    ; mkBond 16 15 17 AromaticBond None
    ; mkBond 17 17 18 SingleBond None
    ; mkBond 18 17 19 AromaticBond None
    ; mkBond 19 11 20 SingleBond None
    ; mkBond 20 20 21 DoubleBond None
    ; mkBond 21 4 22 SingleBond None
    ; mkBond 22 22 23 SingleBond None
    ; mkBond 23 23 1 SingleBond None
    ; mkBond 24 20 8 SingleBond None
    ; mkBond 25 19 12 AromaticBond None
    ; mkBond 26 0 24 SingleBond None
    ; mkBond 27 0 25 SingleBond None
    ; mkBond 28 0 26 SingleBond None
    ; mkBond 29 1 27 SingleBond None
    ; mkBond 30 2 28 SingleBond None
    ; mkBond 31 2 29 SingleBond None
    ; mkBond 32 3 30 SingleBond None
    ; mkBond 33 3 31 SingleBond None
    ; mkBond 34 4 32 SingleBond None
    ; mkBond 35 5 33 SingleBond None
    ; mkBond 36 8 34 SingleBond None
    ; mkBond 37 9 35 SingleBond None
    ; mkBond 38 9 36 SingleBond None
    ; mkBond 39 10 37 SingleBond None
    ; mkBond 40 10 38 SingleBond None
    ; mkBond 41 13 39 SingleBond None
    ; mkBond 42 14 40 SingleBond None
    ; mkBond 43 19 41 SingleBond None
    ; mkBond 44 22 42 SingleBond None
    ; mkBond 45 22 43 SingleBond None
    ; mkBond 46 23 44 SingleBond None
    ; mkBond 47 23 45 SingleBond None
    ].

(* ---- Verification checks ---- *)
(* Confirm demo_molecule is a well-typed Molecule (instant) *)
Check demo_molecule.
(* Confirm property functions are applicable (instant) *)
Check (is_connected demo_molecule).
Check (atom_count demo_molecule).
Check (bond_count demo_molecule).
