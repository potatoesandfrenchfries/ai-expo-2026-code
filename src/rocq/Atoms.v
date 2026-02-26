(** * Atoms: Element types and atomic properties

    Models all 118 chemical elements with their fundamental
    physical and chemical properties as pure Rocq types.
*)

Require Import Stdlib.Reals.Reals.
Require Import Stdlib.ZArith.ZArith.
Require Import Stdlib.Lists.List.
Require Import Stdlib.micromega.Lia.
Require Import Stdlib.micromega.Lra.
Import ListNotations.

Open Scope R_scope.

(** ** 1. Element Types — all 118 elements *)

Inductive Element : Type :=
  (* Period 1 *)
  | eH  | eHe
  (* Period 2 *)
  | eLi | eBe | eB  | eC  | eN  | eO  | eF  | eNe
  (* Period 3 *)
  | eNa | eMg | eAl | eSi | eP  | eS  | eCl | eAr
  (* Period 4 *)
  | eK  | eCa | eSc | eTi | eV  | eCr | eMn | eFe | eCo | eNi | eCu | eZn
  | eGa | eGe | eAs | eSe | eBr | eKr
  (* Period 5 *)
  | eRb | eSr | eY  | eZr | eNb | eMo | eTc | eRu | eRh | ePd | eAg | eCd
  | eIn | eSn | eSb | eTe | eI  | eXe
  (* Period 6 *)
  | eCs | eBa
  (* Lanthanides *)
  | eLa | eCe | ePr | eNd | ePm | eSm | eEu | eGd | eTb | eDy | eHo | eEr | eTm | eYb | eLu
  (* Period 6 d-block *)
  | eHf | eTa | eW  | eRe | eOs | eIr | ePt | eAu | eHg
  | eTl | ePb | eBi | ePo | eAt | eRn
  (* Period 7 *)
  | eFr | eRa
  (* Actinides *)
  | eAc | eTh | ePa | eU  | eNp | ePu | eAm | eCm | eBk | eCf | eEs | eFm | eMd | eNo | eLr
  (* Period 7 d-block *)
  | eRf | eDb | eSg | eBh | eHs | eMt | eDs | eRg | eCn
  | eNh | eFl | eMc | eLv | eTs | eOg.

(** ** 2. Atomic Number *)

Definition atomic_number (e : Element) : nat :=
  match e with
  | eH  => 1  | eHe => 2
  | eLi => 3  | eBe => 4  | eB  => 5  | eC  => 6  | eN  => 7
  | eO  => 8  | eF  => 9  | eNe => 10
  | eNa => 11 | eMg => 12 | eAl => 13 | eSi => 14 | eP  => 15
  | eS  => 16 | eCl => 17 | eAr => 18
  | eK  => 19 | eCa => 20 | eSc => 21 | eTi => 22 | eV  => 23
  | eCr => 24 | eMn => 25 | eFe => 26 | eCo => 27 | eNi => 28
  | eCu => 29 | eZn => 30 | eGa => 31 | eGe => 32 | eAs => 33
  | eSe => 34 | eBr => 35 | eKr => 36
  | eRb => 37 | eSr => 38 | eY  => 39 | eZr => 40 | eNb => 41
  | eMo => 42 | eTc => 43 | eRu => 44 | eRh => 45 | ePd => 46
  | eAg => 47 | eCd => 48 | eIn => 49 | eSn => 50 | eSb => 51
  | eTe => 52 | eI  => 53 | eXe => 54
  | eCs => 55 | eBa => 56
  | eLa => 57 | eCe => 58 | ePr => 59 | eNd => 60 | ePm => 61
  | eSm => 62 | eEu => 63 | eGd => 64 | eTb => 65 | eDy => 66
  | eHo => 67 | eEr => 68 | eTm => 69 | eYb => 70 | eLu => 71
  | eHf => 72 | eTa => 73 | eW  => 74 | eRe => 75 | eOs => 76
  | eIr => 77 | ePt => 78 | eAu => 79 | eHg => 80
  | eTl => 81 | ePb => 82 | eBi => 83 | ePo => 84 | eAt => 85 | eRn => 86
  | eFr => 87 | eRa => 88
  | eAc => 89 | eTh => 90 | ePa => 91 | eU  => 92 | eNp => 93
  | ePu => 94 | eAm => 95 | eCm => 96 | eBk => 97 | eCf => 98
  | eEs => 99 | eFm => 100| eMd => 101| eNo => 102| eLr => 103
  | eRf => 104| eDb => 105| eSg => 106| eBh => 107| eHs => 108
  | eMt => 109| eDs => 110| eRg => 111| eCn => 112
  | eNh => 113| eFl => 114| eMc => 115| eLv => 116| eTs => 117| eOg => 118
  end.

(** ** 3. Atomic Mass (in Daltons, u) *)

Definition atomic_mass (e : Element) : R :=
  match e with
  | eH  => 1.008   | eHe => 4.0026
  | eLi => 6.941   | eBe => 9.0122  | eB  => 10.811  | eC  => 12.011
  | eN  => 14.007  | eO  => 15.999  | eF  => 18.998  | eNe => 20.180
  | eNa => 22.990  | eMg => 24.305  | eAl => 26.982  | eSi => 28.086
  | eP  => 30.974   | eS   => 32.065   | eCl => 35.453   | eAr => 39.948
  | eK  => 39.098  | eCa => 40.078  | eSc => 44.956  | eTi => 47.867
  | eV  => 50.942  | eCr => 51.996  | eMn => 54.938  | eFe => 55.845
  | eCo => 58.933  | eNi => 58.693  | eCu => 63.546  | eZn => 65.38
  | eGa => 69.723  | eGe => 72.630  | eAs => 74.922  | eSe => 78.971
  | eBr => 79.904  | eKr => 83.798
  | eRb => 85.468  | eSr => 87.62   | eY  => 88.906  | eZr => 91.224
  | eNb => 92.906  | eMo => 95.95   | eTc => 98.0    | eRu => 101.07
  | eRh => 102.91  | ePd => 106.42  | eAg => 107.87  | eCd => 112.41
  | eIn => 114.82  | eSn => 118.71  | eSb => 121.76  | eTe => 127.60
  | eI  => 126.90  | eXe => 131.29
  | eCs => 132.91  | eBa => 137.33
  | eLa => 138.91  | eCe => 140.12  | ePr => 140.91  | eNd => 144.24
  | ePm => 145.0   | eSm => 150.36  | eEu => 151.96  | eGd => 157.25
  | eTb => 158.93  | eDy => 162.50  | eHo => 164.93  | eEr => 167.26
  | eTm => 168.93  | eYb => 173.05  | eLu => 174.97
  | eHf => 178.49  | eTa => 180.95  | eW  => 183.84  | eRe => 186.21
  | eOs => 190.23  | eIr => 192.22  | ePt => 195.08  | eAu => 196.97
  | eHg => 200.59
  | eTl => 204.38  | ePb => 207.2   | eBi => 208.98  | ePo => 209.0
  | eAt => 210.0   | eRn => 222.0
  | eFr => 223.0   | eRa => 226.0
  | eAc => 227.0   | eTh => 232.04  | ePa => 231.04  | eU  => 238.03
  | eNp => 237.0   | ePu => 244.0   | eAm => 243.0   | eCm => 247.0
  | eBk => 247.0   | eCf => 251.0   | eEs => 252.0   | eFm => 257.0
  | eMd => 258.0   | eNo => 259.0   | eLr => 262.0
  | eRf => 265.0   | eDb => 268.0   | eSg => 271.0   | eBh => 270.0
  | eHs => 277.0   | eMt => 276.0   | eDs => 281.0   | eRg => 280.0
  | eCn => 285.0   | eNh => 284.0   | eFl => 289.0   | eMc => 288.0
  | eLv => 293.0   | eTs => 294.0   | eOg => 294.0
  end.

(** ** 4. Van der Waals Radius (in Angstroms, Å) *)

Definition van_der_waals_radius (e : Element) : R :=
  match e with
  | eH  => 1.20 | eHe => 1.40
  | eLi => 1.82 | eBe => 1.53 | eB  => 1.92 | eC  => 1.70
  | eN  => 1.55 | eO  => 1.52 | eF  => 1.47 | eNe => 1.54
  | eNa => 2.27 | eMg => 1.73 | eAl => 1.84 | eSi => 2.10
  | eP  => 1.80 | eS  => 1.80 | eCl => 1.75 | eAr => 1.88
  | eK  => 2.75 | eCa => 2.31 | eSc => 2.11 | eTi => 2.00
  | eV  => 2.00 | eCr => 2.00 | eMn => 2.00 | eFe => 2.00
  | eCo => 2.00 | eNi => 1.63 | eCu => 1.40 | eZn => 1.39
  | eGa => 1.87 | eGe => 2.11 | eAs => 1.85 | eSe => 1.90
  | eBr => 1.85 | eKr => 2.02
  | eRb => 3.03 | eSr => 2.49 | eY  => 2.27 | eZr => 2.16
  | eNb => 2.08 | eMo => 2.09 | eTc => 2.16 | eRu => 2.13
  | eRh => 2.10 | ePd => 1.63 | eAg => 1.72 | eCd => 1.58
  | eIn => 1.93 | eSn => 2.17 | eSb => 2.06 | eTe => 2.06
  | eI  => 1.98 | eXe => 2.16
  | eCs => 3.43 | eBa => 2.68
  | eLa => 2.43 | eCe => 2.42 | ePr => 2.40 | eNd => 2.39
  | ePm => 2.38 | eSm => 2.36 | eEu => 2.35 | eGd => 2.34
  | eTb => 2.33 | eDy => 2.31 | eHo => 2.30 | eEr => 2.29
  | eTm => 2.27 | eYb => 2.26 | eLu => 2.24
  | eHf => 2.23 | eTa => 2.22 | eW  => 2.18 | eRe => 2.16
  | eOs => 2.16 | eIr => 2.13 | ePt => 1.75 | eAu => 1.66
  | eHg => 1.55 | eTl => 1.96 | ePb => 2.02 | eBi => 2.07
  | ePo => 1.97 | eAt => 2.02 | eRn => 2.20
  | eFr => 3.48 | eRa => 2.83
  | eAc => 2.47 | eTh => 2.45 | ePa => 2.43 | eU  => 2.41
  | eNp => 2.39 | ePu => 2.43 | eAm => 2.44 | eCm => 2.45
  | eBk => 2.44 | eCf => 2.45 | eEs => 2.45 | eFm => 2.45
  | eMd => 2.46 | eNo => 2.46 | eLr => 2.46
  | eRf => 2.00 | eDb => 2.00 | eSg => 2.00 | eBh => 2.00
  | eHs => 2.00 | eMt => 2.00 | eDs => 2.00 | eRg => 2.00
  | eCn => 2.00 | eNh => 2.00 | eFl => 2.00 | eMc => 2.00
  | eLv => 2.00 | eTs => 2.00 | eOg => 2.00
  end.

(** ** 5. Covalent Radius (in Angstroms, Å) *)

Definition covalent_radius (e : Element) : R :=
  match e with
  | eH  => 0.31 | eHe => 0.28
  | eLi => 1.28 | eBe => 0.96 | eB  => 0.84 | eC  => 0.77
  | eN  => 0.71 | eO  => 0.66 | eF  => 0.57 | eNe => 0.58
  | eNa => 1.66 | eMg => 1.41 | eAl => 1.21 | eSi => 1.11
  | eP  => 1.07 | eS  => 1.05 | eCl => 1.02 | eAr => 1.06
  | eK  => 2.03 | eCa => 1.76 | eSc => 1.70 | eTi => 1.60
  | eV  => 1.53 | eCr => 1.39 | eMn => 1.61 | eFe => 1.52
  | eCo => 1.50 | eNi => 1.24 | eCu => 1.32 | eZn => 1.22
  | eGa => 1.22 | eGe => 1.20 | eAs => 1.19 | eSe => 1.20
  | eBr => 1.20 | eKr => 1.16
  | eRb => 2.20 | eSr => 1.95 | eY  => 1.90 | eZr => 1.75
  | eNb => 1.64 | eMo => 1.54 | eTc => 1.47 | eRu => 1.46
  | eRh => 1.42 | ePd => 1.39 | eAg => 1.45 | eCd => 1.44
  | eIn => 1.42 | eSn => 1.39 | eSb => 1.39 | eTe => 1.38
  | eI  => 1.39 | eXe => 1.40
  | eCs => 2.44 | eBa => 2.15
  | eLa => 2.07 | eCe => 2.04 | ePr => 2.03 | eNd => 2.01
  | ePm => 1.99 | eSm => 1.98 | eEu => 1.98 | eGd => 1.96
  | eTb => 1.94 | eDy => 1.92 | eHo => 1.92 | eEr => 1.89
  | eTm => 1.90 | eYb => 1.87 | eLu => 1.87
  | eHf => 1.75 | eTa => 1.70 | eW  => 1.62 | eRe => 1.51
  | eOs => 1.44 | eIr => 1.41 | ePt => 1.36 | eAu => 1.36
  | eHg => 1.32 | eTl => 1.45 | ePb => 1.46 | eBi => 1.48
  | ePo => 1.40 | eAt => 1.50 | eRn => 1.50
  | eFr => 2.60 | eRa => 2.21
  | eAc => 2.15 | eTh => 2.06 | ePa => 2.00 | eU  => 1.96
  | eNp => 1.90 | ePu => 1.87 | eAm => 1.80 | eCm => 1.69
  | eBk => 1.68 | eCf => 1.68 | eEs => 1.65 | eFm => 1.67
  | eMd => 1.73 | eNo => 1.76 | eLr => 1.61
  | eRf => 1.57 | eDb => 1.49 | eSg => 1.43 | eBh => 1.41
  | eHs => 1.34 | eMt => 1.29 | eDs => 1.28 | eRg => 1.21
  | eCn => 1.22 | eNh => 1.36 | eFl => 1.43 | eMc => 1.62
  | eLv => 1.75 | eTs => 1.65 | eOg => 1.57
  end.

(** ** 6. Pauling Electronegativity *)

(** Electronegativity on the Pauling scale. Noble gases and some heavy
    elements lack well-defined values; we use 0.0 as a sentinel. *)
Definition electronegativity (e : Element) : R :=
  match e with
  | eH  => 2.20 | eHe => 0.00
  | eLi => 0.98 | eBe => 1.57 | eB  => 2.04 | eC  => 2.55
  | eN  => 3.04 | eO  => 3.44 | eF  => 3.98 | eNe => 0.00
  | eNa => 0.93 | eMg => 1.31 | eAl => 1.61 | eSi => 1.90
  | eP  => 2.19 | eS  => 2.58 | eCl => 3.16 | eAr => 0.00
  | eK  => 0.82 | eCa => 1.00 | eSc => 1.36 | eTi => 1.54
  | eV  => 1.63 | eCr => 1.66 | eMn => 1.55 | eFe => 1.83
  | eCo => 1.88 | eNi => 1.91 | eCu => 1.90 | eZn => 1.65
  | eGa => 1.81 | eGe => 2.01 | eAs => 2.18 | eSe => 2.55
  | eBr => 2.96 | eKr => 3.00
  | eRb => 0.82 | eSr => 0.95 | eY  => 1.22 | eZr => 1.33
  | eNb => 1.60 | eMo => 2.16 | eTc => 1.90 | eRu => 2.20
  | eRh => 2.28 | ePd => 2.20 | eAg => 1.93 | eCd => 1.69
  | eIn => 1.78 | eSn => 1.96 | eSb => 2.05 | eTe => 2.10
  | eI  => 2.66 | eXe => 2.60
  | eCs => 0.79 | eBa => 0.89
  | eLa => 1.10 | eCe => 1.12 | ePr => 1.13 | eNd => 1.14
  | ePm => 1.13 | eSm => 1.17 | eEu => 1.20 | eGd => 1.20
  | eTb => 1.10 | eDy => 1.22 | eHo => 1.23 | eEr => 1.24
  | eTm => 1.25 | eYb => 1.10 | eLu => 1.27
  | eHf => 1.30 | eTa => 1.50 | eW  => 2.36 | eRe => 1.90
  | eOs => 2.20 | eIr => 2.20 | ePt => 2.28 | eAu => 2.54
  | eHg => 2.00 | eTl => 1.62 | ePb => 2.33 | eBi => 2.02
  | ePo => 2.00 | eAt => 2.20 | eRn => 0.00
  | eFr => 0.70 | eRa => 0.89
  | eAc => 1.10 | eTh => 1.30 | ePa => 1.50 | eU  => 1.38
  | eNp => 1.36 | ePu => 1.28 | eAm => 1.30 | eCm => 1.30
  | eBk => 1.30 | eCf => 1.30 | eEs => 1.30 | eFm => 1.30
  | eMd => 1.30 | eNo => 1.30 | eLr => 1.30
  | _  => 0.00
  end.

(** ** 7. First Ionization Energy (kJ/mol) *)

Definition ionization_energy (e : Element) : R :=
  match e with
  | eH  => 1312.0 | eHe => 2372.3
  | eLi => 520.2  | eBe => 899.5  | eB  => 800.6  | eC  => 1086.5
  | eN  => 1402.3 | eO  => 1313.9 | eF  => 1681.0 | eNe => 2080.7
  | eNa => 495.8  | eMg => 737.7  | eAl => 577.5  | eSi => 786.5
  | eP  => 1011.8 | eS  => 999.6  | eCl => 1251.2 | eAr => 1520.6
  | eK  => 418.8  | eCa => 589.8  | eSc => 633.1  | eTi => 658.8
  | eV  => 650.9  | eCr => 652.9  | eMn => 717.3  | eFe => 762.5
  | eCo => 760.4  | eNi => 737.1  | eCu => 745.5  | eZn => 906.4
  | eGa => 578.8  | eGe => 762.0  | eAs => 947.0  | eSe => 941.0
  | eBr => 1139.9 | eKr => 1350.8
  | eRb => 403.0  | eSr => 549.5  | eY  => 600.0  | eZr => 640.1
  | eNb => 652.1  | eMo => 684.3  | eTc => 702.0  | eRu => 710.2
  | eRh => 719.7  | ePd => 804.4  | eAg => 731.0  | eCd => 867.8
  | eIn => 558.3  | eSn => 708.6  | eSb => 834.0  | eTe => 869.3
  | eI  => 1008.4 | eXe => 1170.4
  | eCs => 375.7  | eBa => 502.9
  | eLa => 538.1  | eCe => 534.4  | ePr => 527.0  | eNd => 533.1
  | ePm => 540.0  | eSm => 544.5  | eEu => 547.1  | eGd => 593.4
  | eTb => 565.8  | eDy => 573.0  | eHo => 581.0  | eEr => 589.3
  | eTm => 596.7  | eYb => 603.4  | eLu => 523.5
  | eHf => 658.5  | eTa => 761.0  | eW  => 770.0  | eRe => 760.0
  | eOs => 840.0  | eIr => 880.0  | ePt => 870.0  | eAu => 890.1
  | eHg => 1007.1 | eTl => 589.4  | ePb => 715.6  | eBi => 703.0
  | ePo => 812.1  | eAt => 899.0  | eRn => 1037.0
  | eFr => 380.0  | eRa => 509.3
  | eAc => 499.0  | eTh => 587.0  | ePa => 568.0  | eU  => 597.6
  | eNp => 604.5  | ePu => 584.7  | eAm => 578.0  | eCm => 581.0
  | eBk => 601.0  | eCf => 608.0  | eEs => 619.0  | eFm => 629.0
  | eMd => 636.0  | eNo => 641.6  | eLr => 470.0
  | _  => 0.0
  end.

(** ** 8. Electron Affinity (kJ/mol, positive = exothermic) *)

Definition electron_affinity (e : Element) : R :=
  match e with
  | eH  => 72.8   | eHe => 0.0
  | eLi => 59.6   | eBe => 0.0    | eB  => 26.7   | eC  => 121.8
  | eN  => 0.0    | eO  => 141.0  | eF  => 328.2  | eNe => 0.0
  | eNa => 52.9   | eMg => 0.0    | eAl => 41.8   | eSi => 134.1
  | eP  => 72.0   | eS  => 200.4  | eCl => 348.6  | eAr => 0.0
  | eK  => 48.4   | eCa => 2.4    | eSc => 18.1   | eTi => 7.6
  | eV  => 50.7   | eCr => 64.3   | eMn => 0.0    | eFe => 15.7
  | eCo => 63.7   | eNi => 112.0  | eCu => 118.4  | eZn => 0.0
  | eGa => 41.5   | eGe => 119.0  | eAs => 78.2   | eSe => 195.0
  | eBr => 324.6  | eKr => 0.0
  | eRb => 46.9   | eSr => 5.0    | eY  => 29.6   | eZr => 41.1
  | eNb => 86.1   | eMo => 71.9   | eTc => 53.0   | eRu => 101.3
  | eRh => 109.7  | ePd => 53.7   | eAg => 125.6  | eCd => 0.0
  | eIn => 28.9   | eSn => 107.3  | eSb => 101.1  | eTe => 190.2
  | eI  => 295.2  | eXe => 0.0
  | _  => 0.0
  end.

(** ** 9. Period (row in periodic table) *)

Definition period (e : Element) : nat :=
  match e with
  | eH  | eHe => 1
  | eLi | eBe | eB  | eC  | eN  | eO  | eF  | eNe => 2
  | eNa | eMg | eAl | eSi | eP  | eS  | eCl | eAr => 3
  | eK  | eCa | eSc | eTi | eV  | eCr | eMn | eFe | eCo | eNi | eCu | eZn
  | eGa | eGe | eAs | eSe | eBr | eKr => 4
  | eRb | eSr | eY  | eZr | eNb | eMo | eTc | eRu | eRh | ePd | eAg | eCd
  | eIn | eSn | eSb | eTe | eI  | eXe => 5
  | eCs | eBa | eLa | eCe | ePr | eNd | ePm | eSm | eEu | eGd | eTb | eDy
  | eHo | eEr | eTm | eYb | eLu | eHf | eTa | eW  | eRe | eOs | eIr | ePt
  | eAu | eHg | eTl | ePb | eBi | ePo | eAt | eRn => 6
  | eFr | eRa | eAc | eTh | ePa | eU  | eNp | ePu | eAm | eCm | eBk | eCf
  | eEs | eFm | eMd | eNo | eLr | eRf | eDb | eSg | eBh | eHs | eMt | eDs
  | eRg | eCn | eNh | eFl | eMc | eLv | eTs | eOg => 7
  end.

(** ** 10. Group (column in periodic table, using IUPAC numbering 1–18) *)

Definition group (e : Element) : nat :=
  match e with
  | eH  | eLi | eNa | eK  | eRb | eCs | eFr => 1
  | eHe | eNe | eAr | eKr | eXe | eRn | eOg => 18
  | eBe | eMg | eCa | eSr | eBa | eRa => 2
  | eB  | eAl | eGa | eIn | eTl | eNh => 13
  | eC  | eSi | eGe | eSn | ePb | eFl => 14
  | eN  | eP  | eAs | eSb | eBi | eMc => 15
  | eO  | eS  | eSe | eTe | ePo | eLv => 16
  | eF  | eCl | eBr | eI  | eAt | eTs => 17
  | eSc | eY  | eLa | eAc => 3
  | eTi | eZr | eHf | eRf => 4
  | eV  | eNb | eTa | eDb => 5
  | eCr | eMo | eW  | eSg => 6
  | eMn | eTc | eRe | eBh => 7
  | eFe | eRu | eOs | eHs => 8
  | eCo | eRh | eIr | eMt => 9
  | eNi | ePd | ePt | eDs => 10
  | eCu | eAg | eAu | eRg => 11
  | eZn | eCd | eHg | eCn => 12
  (* Lanthanides and actinides: assigned to group 0 as they span 3-17 *)
  | eCe | ePr | eNd | ePm | eSm | eEu | eGd | eTb | eDy | eHo | eEr | eTm | eYb | eLu => 0
  | eTh | ePa | eU  | eNp | ePu | eAm | eCm | eBk | eCf | eEs | eFm | eMd | eNo | eLr => 0
  end.

(** ** 11. Valence Electrons *)

Definition valence_electrons (e : Element) : nat :=
  match e with
  | eH  => 1 | eHe => 2
  | eLi => 1 | eBe => 2 | eB  => 3 | eC  => 4 | eN  => 5 | eO  => 6 | eF  => 7 | eNe => 8
  | eNa => 1 | eMg => 2 | eAl => 3 | eSi => 4 | eP  => 5 | eS  => 6 | eCl => 7 | eAr => 8
  | eK  => 1 | eCa => 2
  | eSc => 3 | eTi => 4 | eV  => 5 | eCr => 6 | eMn => 7 | eFe => 8 | eCo => 9 | eNi => 10
  | eCu => 11| eZn => 12| eGa => 3 | eGe => 4 | eAs => 5 | eSe => 6 | eBr => 7 | eKr => 8
  | eRb => 1 | eSr => 2
  | eY  => 3 | eZr => 4 | eNb => 5 | eMo => 6 | eTc => 7 | eRu => 8 | eRh => 9 | ePd => 10
  | eAg => 11| eCd => 12| eIn => 3 | eSn => 4 | eSb => 5 | eTe => 6 | eI  => 7 | eXe => 8
  | eCs => 1 | eBa => 2
  | eLa => 3 | eCe => 4 | ePr => 5 | eNd => 6 | ePm => 7 | eSm => 8 | eEu => 9 | eGd => 10
  | eTb => 11| eDy => 12| eHo => 13| eEr => 14| eTm => 15| eYb => 16| eLu => 3
  | eHf => 4 | eTa => 5 | eW  => 6 | eRe => 7 | eOs => 8 | eIr => 9 | ePt => 10| eAu => 11
  | eHg => 12| eTl => 3 | ePb => 4 | eBi => 5 | ePo => 6 | eAt => 7 | eRn => 8
  | eFr => 1 | eRa => 2
  | eAc => 3 | eTh => 4 | ePa => 5 | eU  => 6 | eNp => 7 | ePu => 8 | eAm => 9 | eCm => 10
  | eBk => 11| eCf => 12| eEs => 13| eFm => 14| eMd => 15| eNo => 16| eLr => 3
  | eRf => 4 | eDb => 5 | eSg => 6 | eBh => 7 | eHs => 8 | eMt => 9 | eDs => 10| eRg => 11
  | eCn => 12| eNh => 3 | eFl => 4 | eMc => 5 | eLv => 6 | eTs => 7 | eOg => 8
  end.

(** ** 12. Maximum Valence (typical, ignoring hypervalency) *)

Definition max_valence (e : Element) : nat :=
  match e with
  | eH                     => 1
  | eHe | eNe | eAr | eKr     => 0
  | eLi | eNa | eK  | eRb | eCs | eFr => 1
  | eBe | eMg | eCa | eSr | eBa | eRa => 2
  | eB                     => 3
  | eC                     => 4
  | eN                     => 3
  | eO                     => 2
  | eF                     => 1
  | eAl                    => 3
  | eSi                    => 4
  | eP                     => 5   (* hypervalent: PCl₅ *)
  | eS                     => 6   (* hypervalent: SF₆ *)
  | eCl                    => 7   (* ClO₄⁻ *)
  | eBr                    => 5
  | eI                     => 7
  | eXe                    => 8
  | _                     => 6   (* transition metals / heavy elements *)
  end.

(** ** 13. Hypervalent atoms (period ≥ 3 allowing expanded octet) *)

Definition is_hypervalent_possible (e : Element) : bool :=
  match period e with
  | 0 | 1 | 2 => false
  | _          => true
  end.

(** ** 14. Atomic Orbital types present in ground state *)

Inductive OrbitalType : Type := S_orb | P_orb | D_orb | F_orb.

(** Highest occupied orbital type in ground-state configuration *)
Definition highest_orbital (e : Element) : OrbitalType :=
  match e with
  | eH  | eHe => S_orb
  | eLi | eBe => S_orb
  | eB  | eC  | eN  | eO  | eF  | eNe => P_orb
  | eNa | eMg => S_orb
  | eAl | eSi | eP  | eS  | eCl | eAr => P_orb
  | eK  | eCa => S_orb
  | eSc | eTi | eV  | eCr | eMn | eFe | eCo | eNi | eCu | eZn => D_orb
  | eGa | eGe | eAs | eSe | eBr | eKr => P_orb
  | eRb | eSr => S_orb
  | eY  | eZr | eNb | eMo | eTc | eRu | eRh | ePd | eAg | eCd => D_orb
  | eIn | eSn | eSb | eTe | eI  | eXe => P_orb
  | eCs | eBa => S_orb
  | eLa | eCe | ePr | eNd | ePm | eSm | eEu | eGd | eTb | eDy
  | eHo | eEr | eTm | eYb | eLu => F_orb
  | eHf | eTa | eW  | eRe | eOs | eIr | ePt | eAu | eHg => D_orb
  | eTl | ePb | eBi | ePo | eAt | eRn => P_orb
  | eFr | eRa => S_orb
  | eAc | eTh | ePa | eU  | eNp | ePu | eAm | eCm | eBk | eCf
  | eEs | eFm | eMd | eNo | eLr => F_orb
  | eRf | eDb | eSg | eBh | eHs | eMt | eDs | eRg | eCn => D_orb
  | eNh | eFl | eMc | eLv | eTs | eOg => P_orb
  end.

(** ** 15. Element equality decidability *)

Lemma element_eq_dec : forall (a b : Element), {a = b} + {a <> b}.
Proof. decide equality. Defined.

(** ** 16. Basic Properties/Lemmas *)

(** Atomic number is positive for all elements *)
Lemma atomic_number_positive : forall e : Element, (atomic_number e >= 1)%nat.
Proof.
  intro e; destruct e; simpl; lia.
Qed.

(** Atomic number is at most 118 *)
Lemma atomic_number_at_most_118 : forall e : Element, (atomic_number e <= 118)%nat.
Proof.
  intro e; destruct e; simpl; lia.
Qed.

(** Atomic mass is positive *)
Lemma atomic_mass_positive : forall e : Element, atomic_mass e > 0.
Proof.
  intro e; destruct e; simpl; lra.
Qed.

(** Van der Waals radius is positive *)
Lemma vdw_radius_positive : forall e : Element, van_der_waals_radius e > 0.
Proof.
  intro e; destruct e; simpl; lra.
Qed.

(** Electronegativity is non-negative *)
Lemma electronegativity_nonneg : forall e : Element, electronegativity e >= 0.
Proof.
  intro e; destruct e; simpl; lra.
Qed.

(** Fluorine has the highest electronegativity (Pauling scale) *)
Lemma fluorine_max_electronegativity :
  forall e : Element, electronegativity e <= electronegativity eF.
Proof.
  intro e; destruct e; simpl; lra.
Qed.

(** Period is between 1 and 7 *)
Lemma period_bounds : forall e : Element, (1 <= period e /\ period e <= 7)%nat.
Proof.
  intro e; destruct e; simpl; lia.
Qed.

(** Hypervalency requires period >= 3 *)
Lemma hypervalent_requires_period_3 :
  forall e : Element,
    is_hypervalent_possible e = true -> (period e >= 3)%nat.
Proof.
  intros e eH.
  unfold is_hypervalent_possible in eH.
  destruct (period e) eqn:Hp; try discriminate.
  destruct n; try discriminate.
  destruct n; try discriminate.
  lia.
Qed.

Close Scope R_scope.
