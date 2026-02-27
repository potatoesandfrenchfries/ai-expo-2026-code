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

(* Ethanol: CCO with explicit hydrogens — 9 atoms, 8 bonds *)
Definition demo_molecule : Molecule :=
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

(* ---- Verification checks ---- *)
(* Type-check: confirm demo_molecule is a well-typed Molecule *)
Check demo_molecule.

(* Fast runtime evaluation using native OCaml compilation *)
Eval native_compute in (is_connected demo_molecule).
Eval native_compute in (atom_count demo_molecule).
Eval native_compute in (bond_count demo_molecule).
