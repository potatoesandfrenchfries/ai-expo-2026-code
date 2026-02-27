From Coq Require Import List.
Import ListNotations.
Require Import Coq.Reals.Reals.
Open Scope R_scope.

Require Import Coq.ZArith.ZArith.
Require Import Chemistry.Atoms.
Require Import Chemistry.Geometry.
Require Import Chemistry.Bonds.
Require Import Chemistry.Molecules.

Definition demo_molecule : Molecule :=
  mkMol
    [ mkAtom 0 eC (mkPoint (3977 / 1000) (-1442 / 1000) (-33 / 1000)) 0%Z None None ;
      mkAtom 1 eC (mkPoint (3473 / 1000) (-134 / 1000) (544 / 1000)) 0%Z None None ;
      mkAtom 2 eC (mkPoint (4217 / 1000) (1022 / 1000) (-84 / 1000)) 0%Z None None ;
      mkAtom 3 eC (mkPoint (3688 / 1000) (-186 / 1000) (2042 / 1000)) 0%Z None None ;
      mkAtom 4 eC (mkPoint (2011 / 1000) (52 / 1000) (220 / 1000)) 0%Z None None ;
      mkAtom 5 eO (mkPoint (1639 / 1000) (1039 / 1000) (-442 / 1000)) 0%Z None None ;
      mkAtom 6 eN (mkPoint (1044 / 1000) (-858 / 1000) (644 / 1000)) 0%Z None None ;
      mkAtom 7 eC (mkPoint (-376 / 1000) (-808 / 1000) (445 / 1000)) 0%Z None None ;
      mkAtom 8 eS (mkPoint (-1362 / 1000) (-2330 / 1000) (764 / 1000)) 0%Z None None ;
      mkAtom 9 eC (mkPoint (-2969 / 1000) (-1435 / 1000) (302 / 1000)) 0%Z None None ;
      mkAtom 10 eC (mkPoint (-4195 / 1000) (-2281 / 1000) (390 / 1000)) 0%Z None None ;
      mkAtom 11 eC (mkPoint (-4391 / 1000) (-2920 / 1000) (-922 / 1000)) 0%Z None None ;
      mkAtom 12 eN (mkPoint (-5497 / 1000) (-3783 / 1000) (-1103 / 1000)) 0%Z None None ;
      mkAtom 13 eO (mkPoint (-3626 / 1000) (-2738 / 1000) (-1890 / 1000)) 0%Z None None ;
      mkAtom 14 eN (mkPoint (-2565 / 1000) (-426 / 1000) (62 / 1000)) 0%Z None None ;
      mkAtom 15 eC (mkPoint (-1299 / 1000) (48 / 1000) (85 / 1000)) 0%Z None None ;
      mkAtom 16 eC (mkPoint (-1107 / 1000) (1457 / 1000) (-253 / 1000)) 0%Z None None ;
      mkAtom 17 eC (mkPoint (-1429 / 1000) (1983 / 1000) (-1480 / 1000)) 0%Z None None ;
      mkAtom 18 eC (mkPoint (-1342 / 1000) (3329 / 1000) (-1749 / 1000)) 0%Z None None ;
      mkAtom 19 eC (mkPoint (-916 / 1000) (4152 / 1000) (-738 / 1000)) 0%Z None None ;
      mkAtom 20 eC (mkPoint (-576 / 1000) (3689 / 1000) (516 / 1000)) 0%Z None None ;
      mkAtom 21 eF (mkPoint (-161 / 1000) (4557 / 1000) (1470 / 1000)) 0%Z None None ;
      mkAtom 22 eC (mkPoint (-677 / 1000) (2320 / 1000) (751 / 1000)) 0%Z None None ;
      mkAtom 23 eH (mkPoint (3308 / 1000) (-1706 / 1000) (-892 / 1000)) 0%Z None None ;
      mkAtom 24 eH (mkPoint (3782 / 1000) (-2224 / 1000) (742 / 1000)) 0%Z None None ;
      mkAtom 25 eH (mkPoint (5043 / 1000) (-1391 / 1000) (-326 / 1000)) 0%Z None None ;
      mkAtom 26 eH (mkPoint (4959 / 1000) (1445 / 1000) (630 / 1000)) 0%Z None None ;
      mkAtom 27 eH (mkPoint (3519 / 1000) (1829 / 1000) (-420 / 1000)) 0%Z None None ;
      mkAtom 28 eH (mkPoint (4792 / 1000) (688 / 1000) (-972 / 1000)) 0%Z None None ;
      mkAtom 29 eH (mkPoint (2734 / 1000) (-529 / 1000) (2522 / 1000)) 0%Z None None ;
      mkAtom 30 eH (mkPoint (3936 / 1000) (821 / 1000) (2437 / 1000)) 0%Z None None ;
      mkAtom 31 eH (mkPoint (4469 / 1000) (-931 / 1000) (2280 / 1000)) 0%Z None None ;
      mkAtom 32 eH (mkPoint (1450 / 1000) (-1691 / 1000) (1193 / 1000)) 0%Z None None ;
      mkAtom 33 eH (mkPoint (-5017 / 1000) (-1613 / 1000) (688 / 1000)) 0%Z None None ;
      mkAtom 34 eH (mkPoint (-4113 / 1000) (-3053 / 1000) (1189 / 1000)) 0%Z None None ;
      mkAtom 35 eH (mkPoint (-5366 / 1000) (-4800 / 1000) (-1216 / 1000)) 0%Z None None ;
      mkAtom 36 eH (mkPoint (-6449 / 1000) (-3371 / 1000) (-1120 / 1000)) 0%Z None None ;
      mkAtom 37 eH (mkPoint (-1760 / 1000) (1298 / 1000) (-2247 / 1000)) 0%Z None None ;
      mkAtom 38 eH (mkPoint (-1596 / 1000) (3734 / 1000) (-2710 / 1000)) 0%Z None None ;
      mkAtom 39 eH (mkPoint (-840 / 1000) (5208 / 1000) (-927 / 1000)) 0%Z None None ;
      mkAtom 40 eH (mkPoint (-407 / 1000) (1978 / 1000) (1738 / 1000)) 0%Z None None ]
    [].
