From Stdlib Require Import List.
Import ListNotations.
Require Import Stdlib.Reals.Reals.
Open Scope R_scope.

Require Import Stdlib.ZArith.ZArith.
Require Import Stdlib.Arith.Arith.
Require Import Stdlib.micromega.Lia.
Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.
Require Import Chemistry.Molecules.

Definition demo_molecule : Molecule :=
  mkMol
    [ mkAtom 0 eO (mkPoint (1807 / 1000) (1095 / 1000) (-1668 / 1000)) 0%Z None None ;
      mkAtom 1 eC (mkPoint (1243 / 1000) (561 / 1000) (-688 / 1000)) 0%Z None None ;
      mkAtom 2 eN (mkPoint (1997 / 1000) (-67 / 1000) (294 / 1000)) 0%Z None None ;
      mkAtom 3 eC (mkPoint (3411 / 1000) (-201 / 1000) (353 / 1000)) 0%Z None None ;
      mkAtom 4 eC (mkPoint (3990 / 1000) (-877 / 1000) (1437 / 1000)) 0%Z None None ;
      mkAtom 5 eC (mkPoint (5346 / 1000) (-1043 / 1000) (1558 / 1000)) 0%Z None None ;
      mkAtom 6 eC (mkPoint (6196 / 1000) (-530 / 1000) (582 / 1000)) 0%Z None None ;
      mkAtom 7 eC (mkPoint (5668 / 1000) (133 / 1000) (-486 / 1000)) 0%Z None None ;
      mkAtom 8 eCl (mkPoint (6648 / 1000) (815 / 1000) (-1765 / 1000)) 0%Z None None ;
      mkAtom 9 eC (mkPoint (4296 / 1000) (278 / 1000) (-570 / 1000)) 0%Z None None ;
      mkAtom 10 eC (mkPoint (-218 / 1000) (579 / 1000) (-543 / 1000)) 0%Z None None ;
      mkAtom 11 eS (mkPoint (-1315 / 1000) (1346 / 1000) (-1721 / 1000)) 0%Z None None ;
      mkAtom 12 eC (mkPoint (-2774 / 1000) (986 / 1000) (-932 / 1000)) 0%Z None None ;
      mkAtom 13 eC (mkPoint (-2391 / 1000) (302 / 1000) (189 / 1000)) 0%Z None None ;
      mkAtom 14 eC (mkPoint (-1068 / 1000) (106 / 1000) (366 / 1000)) 0%Z None None ;
      mkAtom 15 eC (mkPoint (-3570 / 1000) (-121 / 1000) (1055 / 1000)) 0%Z None None ;
      mkAtom 16 eC (mkPoint (-4607 / 1000) (-811 / 1000) (208 / 1000)) 0%Z None None ;
      mkAtom 17 eC (mkPoint (-4890 / 1000) (-133 / 1000) (-1085 / 1000)) 0%Z None None ;
      mkAtom 18 eC (mkPoint (-4254 / 1000) (1239 / 1000) (-1186 / 1000)) 0%Z None None ;
      mkAtom 19 eH (mkPoint (1422 / 1000) (-496 / 1000) (1096 / 1000)) 0%Z None None ;
      mkAtom 20 eH (mkPoint (3333 / 1000) (-1284 / 1000) (2209 / 1000)) 0%Z None None ;
      mkAtom 21 eH (mkPoint (5788 / 1000) (-1565 / 1000) (2395 / 1000)) 0%Z None None ;
      mkAtom 22 eH (mkPoint (7275 / 1000) (-641 / 1000) (642 / 1000)) 0%Z None None ;
      mkAtom 23 eH (mkPoint (3938 / 1000) (805 / 1000) (-1425 / 1000)) 0%Z None None ;
      mkAtom 24 eH (mkPoint (-670 / 1000) (-444 / 1000) (1259 / 1000)) 0%Z None None ;
      mkAtom 25 eH (mkPoint (-3247 / 1000) (-851 / 1000) (1825 / 1000)) 0%Z None None ;
      mkAtom 26 eH (mkPoint (-3954 / 1000) (801 / 1000) (1536 / 1000)) 0%Z None None ;
      mkAtom 27 eH (mkPoint (-4362 / 1000) (-1875 / 1000) (76 / 1000)) 0%Z None None ;
      mkAtom 28 eH (mkPoint (-5560 / 1000) (-796 / 1000) (813 / 1000)) 0%Z None None ;
      mkAtom 29 eH (mkPoint (-5986 / 1000) (-60 / 1000) (-1292 / 1000)) 0%Z None None ;
      mkAtom 30 eH (mkPoint (-4453 / 1000) (-766 / 1000) (-1901 / 1000)) 0%Z None None ;
      mkAtom 31 eH (mkPoint (-4673 / 1000) (1911 / 1000) (-436 / 1000)) 0%Z None None ;
      mkAtom 32 eH (mkPoint (-4363 / 1000) (1607 / 1000) (-2224 / 1000)) 0%Z None None ]
    [ mkBond 0 0 1 DoubleBond None ;
      mkBond 1 1 2 SingleBond None ;
      mkBond 2 2 3 SingleBond None ;
      mkBond 3 3 4 AromaticBond None ;
      mkBond 4 4 5 AromaticBond None ;
      mkBond 5 5 6 AromaticBond None ;
      mkBond 6 6 7 AromaticBond None ;
      mkBond 7 7 8 SingleBond None ;
      mkBond 8 7 9 AromaticBond None ;
      mkBond 9 1 10 SingleBond None ;
      mkBond 10 10 11 AromaticBond None ;
      mkBond 11 11 12 AromaticBond None ;
      mkBond 12 12 13 AromaticBond None ;
      mkBond 13 13 14 AromaticBond None ;
      mkBond 14 13 15 SingleBond None ;
      mkBond 15 15 16 SingleBond None ;
      mkBond 16 16 17 SingleBond None ;
      mkBond 17 17 18 SingleBond None ;
      mkBond 18 9 3 AromaticBond None ;
      mkBond 19 14 10 AromaticBond None ;
      mkBond 20 18 12 SingleBond None ;
      mkBond 21 2 19 SingleBond None ;
      mkBond 22 4 20 SingleBond None ;
      mkBond 23 5 21 SingleBond None ;
      mkBond 24 6 22 SingleBond None ;
      mkBond 25 9 23 SingleBond None ;
      mkBond 26 14 24 SingleBond None ;
      mkBond 27 15 25 SingleBond None ;
      mkBond 28 15 26 SingleBond None ;
      mkBond 29 16 27 SingleBond None ;
      mkBond 30 16 28 SingleBond None ;
      mkBond 31 17 29 SingleBond None ;
      mkBond 32 17 30 SingleBond None ;
      mkBond 33 18 31 SingleBond None ;
      mkBond 34 18 32 SingleBond None ].

Theorem generated_has_atoms : length (mol_atoms demo_molecule) > 0.
Proof. compute; lia. Qed.

Theorem generated_has_bonds : length (mol_bonds demo_molecule) > 0.
Proof. compute; lia. Qed.
