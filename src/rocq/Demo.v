From Stdlib Require Import List.
Import ListNotations.
Require Import Stdlib.Reals.Reals.
Open Scope R_scope.

Require Import Stdlib.ZArith.ZArith.
Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.
Require Import Chemistry.Molecules.

Definition demo_molecule : Molecule :=
  mkMol
    [ mkAtom 0 eC (mkPoint (6055 / 1000) (1693 / 1000) (198 / 1000)) 0%Z None None ;
      mkAtom 1 eO (mkPoint (4839 / 1000) (2168 / 1000) (736 / 1000)) 0%Z None None ;
      mkAtom 2 eC (mkPoint (3742 / 1000) (1346 / 1000) (786 / 1000)) 0%Z None None ;
      mkAtom 3 eC (mkPoint (3824 / 1000) (58 / 1000) (313 / 1000)) 0%Z None None ;
      mkAtom 4 eC (mkPoint (2722 / 1000) (-787 / 1000) (356 / 1000)) 0%Z None None ;
      mkAtom 5 eC (mkPoint (1525 / 1000) (-333 / 1000) (880 / 1000)) 0%Z None None ;
      mkAtom 6 eS (mkPoint (121 / 1000) (-1400 / 1000) (939 / 1000)) 0%Z None None ;
      mkAtom 7 eO (mkPoint (-733 / 1000) (-1120 / 1000) (2124 / 1000)) 0%Z None None ;
      mkAtom 8 eO (mkPoint (578 / 1000) (-2845 / 1000) (1003 / 1000)) 0%Z None None ;
      mkAtom 9 eN (mkPoint (-752 / 1000) (-1208 / 1000) (-517 / 1000)) 0%Z None None ;
      mkAtom 10 eC (mkPoint (-1224 / 1000) (146 / 1000) (-665 / 1000)) 0%Z None None ;
      mkAtom 11 eC (mkPoint (-2375 / 1000) (125 / 1000) (-1643 / 1000)) 0%Z None None ;
      mkAtom 12 eC (mkPoint (-3581 / 1000) (-449 / 1000) (-920 / 1000)) 0%Z None None ;
      mkAtom 13 eC (mkPoint (-4086 / 1000) (628 / 1000) (-19 / 1000)) 0%Z None None ;
      mkAtom 14 eN (mkPoint (-4493 / 1000) (1835 / 1000) (-654 / 1000)) 0%Z None None ;
      mkAtom 15 eO (mkPoint (-4146 / 1000) (476 / 1000) (1222 / 1000)) 0%Z None None ;
      mkAtom 16 eC (mkPoint (-3137 / 1000) (-1658 / 1000) (-157 / 1000)) 0%Z None None ;
      mkAtom 17 eC (mkPoint (-1793 / 1000) (-2187 / 1000) (-619 / 1000)) 0%Z None None ;
      mkAtom 18 eC (mkPoint (1474 / 1000) (963 / 1000) (1347 / 1000)) 0%Z None None ;
      mkAtom 19 eC (mkPoint (2561 / 1000) (1810 / 1000) (1308 / 1000)) 0%Z None None ;
      mkAtom 20 eH (mkPoint (6483 / 1000) (856 / 1000) (817 / 1000)) 0%Z None None ;
      mkAtom 21 eH (mkPoint (6797 / 1000) (2510 / 1000) (216 / 1000)) 0%Z None None ;
      mkAtom 22 eH (mkPoint (5920 / 1000) (1279 / 1000) (-812 / 1000)) 0%Z None None ;
      mkAtom 23 eH (mkPoint (4755 / 1000) (-300 / 1000) (-96 / 1000)) 0%Z None None ;
      mkAtom 24 eH (mkPoint (2777 / 1000) (-1809 / 1000) (-15 / 1000)) 0%Z None None ;
      mkAtom 25 eH (mkPoint (-1575 / 1000) (500 / 1000) (332 / 1000)) 0%Z None None ;
      mkAtom 26 eH (mkPoint (-413 / 1000) (821 / 1000) (-1013 / 1000)) 0%Z None None ;
      mkAtom 27 eH (mkPoint (-2668 / 1000) (1164 / 1000) (-1942 / 1000)) 0%Z None None ;
      mkAtom 28 eH (mkPoint (-2080 / 1000) (-488 / 1000) (-2506 / 1000)) 0%Z None None ;
      mkAtom 29 eH (mkPoint (-4407 / 1000) (-673 / 1000) (-1630 / 1000)) 0%Z None None ;
      mkAtom 30 eH (mkPoint (-4478 / 1000) (2764 / 1000) (-159 / 1000)) 0%Z None None ;
      mkAtom 31 eH (mkPoint (-4825 / 1000) (1852 / 1000) (-1642 / 1000)) 0%Z None None ;
      mkAtom 32 eH (mkPoint (-3879 / 1000) (-2478 / 1000) (-371 / 1000)) 0%Z None None ;
      mkAtom 33 eH (mkPoint (-3086 / 1000) (-1507 / 1000) (933 / 1000)) 0%Z None None ;
      mkAtom 34 eH (mkPoint (-1877 / 1000) (-2642 / 1000) (-1630 / 1000)) 0%Z None None ;
      mkAtom 35 eH (mkPoint (-1561 / 1000) (-3051 / 1000) (67 / 1000)) 0%Z None None ;
      mkAtom 36 eH (mkPoint (511 / 1000) (1290 / 1000) (1755 / 1000)) 0%Z None None ;
      mkAtom 37 eH (mkPoint (2486 / 1000) (2816 / 1000) (1681 / 1000)) 0%Z None None ]
    [].
