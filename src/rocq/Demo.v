From Stdlib Require Import List.
Import ListNotations.
Require Import Stdlib.Reals.Reals.
Open Scope R_scope.

Require Import Stdlib.ZArith.ZArith.
Require Import Stdlib.ZArith.ZArith.
Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.
Require Import Chemistry.Molecule.

Definition demo_molecule : Molecule :=
  mkMolecule
    [ mkAtomInstance eC (mkPoint (-6405 / 1000) (914 / 1000) (512 / 1000)) ;
      mkAtomInstance eN (mkPoint (-5546 / 1000) (-164 / 1000) (42 / 1000)) ;
      mkAtomInstance eC (mkPoint (-4674 / 1000) (-500 / 1000) (1145 / 1000)) ]
    [].
