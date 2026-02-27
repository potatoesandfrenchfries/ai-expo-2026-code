(** * Molecule: Definition of a Chemical Molecule *)

From Stdlib Require Import List.
Import ListNotations.

Require Import Chemistry.Atoms.
Require Import Chemistry.Bonds.
Require Import Chemistry.Geometry.

(** An atom instance in a molecule has an element type and a 3D position. *)
Record AtomInstance : Type := mkAtomInstance
  { element : Element ;
    position : Point3D }.

(** A bond connects two atoms (by their indices in the atom list) with a specific bond type. *)
Record BondInstance : Type := mkBondInstance
  { atom1_idx : nat ;
    atom2_idx : nat ;
    btype : BondType }.

(** A molecule is a collection of atoms and the bonds between them. *)
Record Molecule : Type := mkMolecule
  { atoms : list AtomInstance ;
    bonds : list BondInstance }.
