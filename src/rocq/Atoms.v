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
  | H  | He
  (* Period 2 *)
  | Li | Be | B  | C  | N  | O  | F  | Ne
  (* Period 3 *)
  | Na | Mg | Al | Si | P  | S  | Cl | Ar
  (* Period 4 *)
  | K  | Ca | Sc | Ti | V  | Cr | Mn | Fe | Co | Ni | Cu | Zn
  | Ga | Ge | As | Se | Br | Kr
  (* Period 5 *)
  | Rb | Sr | Y  | Zr | Nb | Mo | Tc | Ru | Rh | Pd | Ag | Cd
  | In | Sn | Sb | Te | I  | Xe
  (* Period 6 *)
  | Cs | Ba
  (* Lanthanides *)
  | La | Ce | Pr | Nd | Pm | Sm | Eu | Gd | Tb | Dy | Ho | Er | Tm | Yb | Lu
  (* Period 6 d-block *)
  | Hf | Ta | W  | Re | Os | Ir | Pt | Au | Hg
  | Tl | Pb | Bi | Po | At | Rn
  (* Period 7 *)
  | Fr | Ra
  (* Actinides *)
  | Ac | Th | Pa | U  | Np | Pu | Am | Cm | Bk | Cf | Es | Fm | Md | No | Lr
  (* Period 7 d-block *)
  | Rf | Db | Sg | Bh | Hs | Mt | Ds | Rg | Cn
  | Nh | Fl | Mc | Lv | Ts | Og.

(** ** 2. Atomic Number *)

Definition atomic_number (e : Element) : nat :=
  match e with
  | H  => 1  | He => 2
  | Li => 3  | Be => 4  | B  => 5  | C  => 6  | N  => 7
  | O  => 8  | F  => 9  | Ne => 10
  | Na => 11 | Mg => 12 | Al => 13 | Si => 14 | P  => 15
  | S  => 16 | Cl => 17 | Ar => 18
  | K  => 19 | Ca => 20 | Sc => 21 | Ti => 22 | V  => 23
  | Cr => 24 | Mn => 25 | Fe => 26 | Co => 27 | Ni => 28
  | Cu => 29 | Zn => 30 | Ga => 31 | Ge => 32 | As => 33
  | Se => 34 | Br => 35 | Kr => 36
  | Rb => 37 | Sr => 38 | Y  => 39 | Zr => 40 | Nb => 41
  | Mo => 42 | Tc => 43 | Ru => 44 | Rh => 45 | Pd => 46
  | Ag => 47 | Cd => 48 | In => 49 | Sn => 50 | Sb => 51
  | Te => 52 | I  => 53 | Xe => 54
  | Cs => 55 | Ba => 56
  | La => 57 | Ce => 58 | Pr => 59 | Nd => 60 | Pm => 61
  | Sm => 62 | Eu => 63 | Gd => 64 | Tb => 65 | Dy => 66
  | Ho => 67 | Er => 68 | Tm => 69 | Yb => 70 | Lu => 71
  | Hf => 72 | Ta => 73 | W  => 74 | Re => 75 | Os => 76
  | Ir => 77 | Pt => 78 | Au => 79 | Hg => 80
  | Tl => 81 | Pb => 82 | Bi => 83 | Po => 84 | At => 85 | Rn => 86
  | Fr => 87 | Ra => 88
  | Ac => 89 | Th => 90 | Pa => 91 | U  => 92 | Np => 93
  | Pu => 94 | Am => 95 | Cm => 96 | Bk => 97 | Cf => 98
  | Es => 99 | Fm => 100| Md => 101| No => 102| Lr => 103
  | Rf => 104| Db => 105| Sg => 106| Bh => 107| Hs => 108
  | Mt => 109| Ds => 110| Rg => 111| Cn => 112
  | Nh => 113| Fl => 114| Mc => 115| Lv => 116| Ts => 117| Og => 118
  end.

(** ** 3. Atomic Mass (in Daltons, u) *)

Definition atomic_mass (e : Element) : R :=
  match e with
  | H  => 1.008   | He => 4.0026
  | Li => 6.941   | Be => 9.0122  | B  => 10.811  | C  => 12.011
  | N  => 14.007  | O  => 15.999  | F  => 18.998  | Ne => 20.180
  | Na => 22.990  | Mg => 24.305  | Al => 26.982  | Si => 28.086
  | P  => 30.974  | S  => 32.065  | Cl => 35.453  | Ar => 39.948
  | K  => 39.098  | Ca => 40.078  | Sc => 44.956  | Ti => 47.867
  | V  => 50.942  | Cr => 51.996  | Mn => 54.938  | Fe => 55.845
  | Co => 58.933  | Ni => 58.693  | Cu => 63.546  | Zn => 65.38
  | Ga => 69.723  | Ge => 72.630  | As => 74.922  | Se => 78.971
  | Br => 79.904  | Kr => 83.798
  | Rb => 85.468  | Sr => 87.62   | Y  => 88.906  | Zr => 91.224
  | Nb => 92.906  | Mo => 95.95   | Tc => 98.0    | Ru => 101.07
  | Rh => 102.91  | Pd => 106.42  | Ag => 107.87  | Cd => 112.41
  | In => 114.82  | Sn => 118.71  | Sb => 121.76  | Te => 127.60
  | I  => 126.90  | Xe => 131.29
  | Cs => 132.91  | Ba => 137.33
  | La => 138.91  | Ce => 140.12  | Pr => 140.91  | Nd => 144.24
  | Pm => 145.0   | Sm => 150.36  | Eu => 151.96  | Gd => 157.25
  | Tb => 158.93  | Dy => 162.50  | Ho => 164.93  | Er => 167.26
  | Tm => 168.93  | Yb => 173.05  | Lu => 174.97
  | Hf => 178.49  | Ta => 180.95  | W  => 183.84  | Re => 186.21
  | Os => 190.23  | Ir => 192.22  | Pt => 195.08  | Au => 196.97
  | Hg => 200.59
  | Tl => 204.38  | Pb => 207.2   | Bi => 208.98  | Po => 209.0
  | At => 210.0   | Rn => 222.0
  | Fr => 223.0   | Ra => 226.0
  | Ac => 227.0   | Th => 232.04  | Pa => 231.04  | U  => 238.03
  | Np => 237.0   | Pu => 244.0   | Am => 243.0   | Cm => 247.0
  | Bk => 247.0   | Cf => 251.0   | Es => 252.0   | Fm => 257.0
  | Md => 258.0   | No => 259.0   | Lr => 262.0
  | Rf => 265.0   | Db => 268.0   | Sg => 271.0   | Bh => 270.0
  | Hs => 277.0   | Mt => 276.0   | Ds => 281.0   | Rg => 280.0
  | Cn => 285.0   | Nh => 284.0   | Fl => 289.0   | Mc => 288.0
  | Lv => 293.0   | Ts => 294.0   | Og => 294.0
  end.

(** ** 4. Van der Waals Radius (in Angstroms, Å) *)

Definition van_der_waals_radius (e : Element) : R :=
  match e with
  | H  => 1.20 | He => 1.40
  | Li => 1.82 | Be => 1.53 | B  => 1.92 | C  => 1.70
  | N  => 1.55 | O  => 1.52 | F  => 1.47 | Ne => 1.54
  | Na => 2.27 | Mg => 1.73 | Al => 1.84 | Si => 2.10
  | P  => 1.80 | S  => 1.80 | Cl => 1.75 | Ar => 1.88
  | K  => 2.75 | Ca => 2.31 | Sc => 2.11 | Ti => 2.00
  | V  => 2.00 | Cr => 2.00 | Mn => 2.00 | Fe => 2.00
  | Co => 2.00 | Ni => 1.63 | Cu => 1.40 | Zn => 1.39
  | Ga => 1.87 | Ge => 2.11 | As => 1.85 | Se => 1.90
  | Br => 1.85 | Kr => 2.02
  | Rb => 3.03 | Sr => 2.49 | Y  => 2.27 | Zr => 2.16
  | Nb => 2.08 | Mo => 2.09 | Tc => 2.16 | Ru => 2.13
  | Rh => 2.10 | Pd => 1.63 | Ag => 1.72 | Cd => 1.58
  | In => 1.93 | Sn => 2.17 | Sb => 2.06 | Te => 2.06
  | I  => 1.98 | Xe => 2.16
  | Cs => 3.43 | Ba => 2.68
  | La => 2.43 | Ce => 2.42 | Pr => 2.40 | Nd => 2.39
  | Pm => 2.38 | Sm => 2.36 | Eu => 2.35 | Gd => 2.34
  | Tb => 2.33 | Dy => 2.31 | Ho => 2.30 | Er => 2.29
  | Tm => 2.27 | Yb => 2.26 | Lu => 2.24
  | Hf => 2.23 | Ta => 2.22 | W  => 2.18 | Re => 2.16
  | Os => 2.16 | Ir => 2.13 | Pt => 1.75 | Au => 1.66
  | Hg => 1.55 | Tl => 1.96 | Pb => 2.02 | Bi => 2.07
  | Po => 1.97 | At => 2.02 | Rn => 2.20
  | Fr => 3.48 | Ra => 2.83
  | Ac => 2.47 | Th => 2.45 | Pa => 2.43 | U  => 2.41
  | Np => 2.39 | Pu => 2.43 | Am => 2.44 | Cm => 2.45
  | Bk => 2.44 | Cf => 2.45 | Es => 2.45 | Fm => 2.45
  | Md => 2.46 | No => 2.46 | Lr => 2.46
  | Rf => 2.00 | Db => 2.00 | Sg => 2.00 | Bh => 2.00
  | Hs => 2.00 | Mt => 2.00 | Ds => 2.00 | Rg => 2.00
  | Cn => 2.00 | Nh => 2.00 | Fl => 2.00 | Mc => 2.00
  | Lv => 2.00 | Ts => 2.00 | Og => 2.00
  end.

(** ** 5. Covalent Radius (in Angstroms, Å) *)

Definition covalent_radius (e : Element) : R :=
  match e with
  | H  => 0.31 | He => 0.28
  | Li => 1.28 | Be => 0.96 | B  => 0.84 | C  => 0.77
  | N  => 0.71 | O  => 0.66 | F  => 0.57 | Ne => 0.58
  | Na => 1.66 | Mg => 1.41 | Al => 1.21 | Si => 1.11
  | P  => 1.07 | S  => 1.05 | Cl => 1.02 | Ar => 1.06
  | K  => 2.03 | Ca => 1.76 | Sc => 1.70 | Ti => 1.60
  | V  => 1.53 | Cr => 1.39 | Mn => 1.61 | Fe => 1.52
  | Co => 1.50 | Ni => 1.24 | Cu => 1.32 | Zn => 1.22
  | Ga => 1.22 | Ge => 1.20 | As => 1.19 | Se => 1.20
  | Br => 1.20 | Kr => 1.16
  | Rb => 2.20 | Sr => 1.95 | Y  => 1.90 | Zr => 1.75
  | Nb => 1.64 | Mo => 1.54 | Tc => 1.47 | Ru => 1.46
  | Rh => 1.42 | Pd => 1.39 | Ag => 1.45 | Cd => 1.44
  | In => 1.42 | Sn => 1.39 | Sb => 1.39 | Te => 1.38
  | I  => 1.39 | Xe => 1.40
  | Cs => 2.44 | Ba => 2.15
  | La => 2.07 | Ce => 2.04 | Pr => 2.03 | Nd => 2.01
  | Pm => 1.99 | Sm => 1.98 | Eu => 1.98 | Gd => 1.96
  | Tb => 1.94 | Dy => 1.92 | Ho => 1.92 | Er => 1.89
  | Tm => 1.90 | Yb => 1.87 | Lu => 1.87
  | Hf => 1.75 | Ta => 1.70 | W  => 1.62 | Re => 1.51
  | Os => 1.44 | Ir => 1.41 | Pt => 1.36 | Au => 1.36
  | Hg => 1.32 | Tl => 1.45 | Pb => 1.46 | Bi => 1.48
  | Po => 1.40 | At => 1.50 | Rn => 1.50
  | Fr => 2.60 | Ra => 2.21
  | Ac => 2.15 | Th => 2.06 | Pa => 2.00 | U  => 1.96
  | Np => 1.90 | Pu => 1.87 | Am => 1.80 | Cm => 1.69
  | Bk => 1.68 | Cf => 1.68 | Es => 1.65 | Fm => 1.67
  | Md => 1.73 | No => 1.76 | Lr => 1.61
  | Rf => 1.57 | Db => 1.49 | Sg => 1.43 | Bh => 1.41
  | Hs => 1.34 | Mt => 1.29 | Ds => 1.28 | Rg => 1.21
  | Cn => 1.22 | Nh => 1.36 | Fl => 1.43 | Mc => 1.62
  | Lv => 1.75 | Ts => 1.65 | Og => 1.57
  end.

(** ** 6. Pauling Electronegativity *)

(** Electronegativity on the Pauling scale. Noble gases and some heavy
    elements lack well-defined values; we use 0.0 as a sentinel. *)
Definition electronegativity (e : Element) : R :=
  match e with
  | H  => 2.20 | He => 0.00
  | Li => 0.98 | Be => 1.57 | B  => 2.04 | C  => 2.55
  | N  => 3.04 | O  => 3.44 | F  => 3.98 | Ne => 0.00
  | Na => 0.93 | Mg => 1.31 | Al => 1.61 | Si => 1.90
  | P  => 2.19 | S  => 2.58 | Cl => 3.16 | Ar => 0.00
  | K  => 0.82 | Ca => 1.00 | Sc => 1.36 | Ti => 1.54
  | V  => 1.63 | Cr => 1.66 | Mn => 1.55 | Fe => 1.83
  | Co => 1.88 | Ni => 1.91 | Cu => 1.90 | Zn => 1.65
  | Ga => 1.81 | Ge => 2.01 | As => 2.18 | Se => 2.55
  | Br => 2.96 | Kr => 3.00
  | Rb => 0.82 | Sr => 0.95 | Y  => 1.22 | Zr => 1.33
  | Nb => 1.60 | Mo => 2.16 | Tc => 1.90 | Ru => 2.20
  | Rh => 2.28 | Pd => 2.20 | Ag => 1.93 | Cd => 1.69
  | In => 1.78 | Sn => 1.96 | Sb => 2.05 | Te => 2.10
  | I  => 2.66 | Xe => 2.60
  | Cs => 0.79 | Ba => 0.89
  | La => 1.10 | Ce => 1.12 | Pr => 1.13 | Nd => 1.14
  | Pm => 1.13 | Sm => 1.17 | Eu => 1.20 | Gd => 1.20
  | Tb => 1.10 | Dy => 1.22 | Ho => 1.23 | Er => 1.24
  | Tm => 1.25 | Yb => 1.10 | Lu => 1.27
  | Hf => 1.30 | Ta => 1.50 | W  => 2.36 | Re => 1.90
  | Os => 2.20 | Ir => 2.20 | Pt => 2.28 | Au => 2.54
  | Hg => 2.00 | Tl => 1.62 | Pb => 2.33 | Bi => 2.02
  | Po => 2.00 | At => 2.20 | Rn => 0.00
  | Fr => 0.70 | Ra => 0.89
  | Ac => 1.10 | Th => 1.30 | Pa => 1.50 | U  => 1.38
  | Np => 1.36 | Pu => 1.28 | Am => 1.30 | Cm => 1.30
  | Bk => 1.30 | Cf => 1.30 | Es => 1.30 | Fm => 1.30
  | Md => 1.30 | No => 1.30 | Lr => 1.30
  | _  => 0.00
  end.

(** ** 7. First Ionization Energy (kJ/mol) *)

Definition ionization_energy (e : Element) : R :=
  match e with
  | H  => 1312.0 | He => 2372.3
  | Li => 520.2  | Be => 899.5  | B  => 800.6  | C  => 1086.5
  | N  => 1402.3 | O  => 1313.9 | F  => 1681.0 | Ne => 2080.7
  | Na => 495.8  | Mg => 737.7  | Al => 577.5  | Si => 786.5
  | P  => 1011.8 | S  => 999.6  | Cl => 1251.2 | Ar => 1520.6
  | K  => 418.8  | Ca => 589.8  | Sc => 633.1  | Ti => 658.8
  | V  => 650.9  | Cr => 652.9  | Mn => 717.3  | Fe => 762.5
  | Co => 760.4  | Ni => 737.1  | Cu => 745.5  | Zn => 906.4
  | Ga => 578.8  | Ge => 762.0  | As => 947.0  | Se => 941.0
  | Br => 1139.9 | Kr => 1350.8
  | Rb => 403.0  | Sr => 549.5  | Y  => 600.0  | Zr => 640.1
  | Nb => 652.1  | Mo => 684.3  | Tc => 702.0  | Ru => 710.2
  | Rh => 719.7  | Pd => 804.4  | Ag => 731.0  | Cd => 867.8
  | In => 558.3  | Sn => 708.6  | Sb => 834.0  | Te => 869.3
  | I  => 1008.4 | Xe => 1170.4
  | Cs => 375.7  | Ba => 502.9
  | La => 538.1  | Ce => 534.4  | Pr => 527.0  | Nd => 533.1
  | Pm => 540.0  | Sm => 544.5  | Eu => 547.1  | Gd => 593.4
  | Tb => 565.8  | Dy => 573.0  | Ho => 581.0  | Er => 589.3
  | Tm => 596.7  | Yb => 603.4  | Lu => 523.5
  | Hf => 658.5  | Ta => 761.0  | W  => 770.0  | Re => 760.0
  | Os => 840.0  | Ir => 880.0  | Pt => 870.0  | Au => 890.1
  | Hg => 1007.1 | Tl => 589.4  | Pb => 715.6  | Bi => 703.0
  | Po => 812.1  | At => 899.0  | Rn => 1037.0
  | Fr => 380.0  | Ra => 509.3
  | Ac => 499.0  | Th => 587.0  | Pa => 568.0  | U  => 597.6
  | Np => 604.5  | Pu => 584.7  | Am => 578.0  | Cm => 581.0
  | Bk => 601.0  | Cf => 608.0  | Es => 619.0  | Fm => 629.0
  | Md => 636.0  | No => 641.6  | Lr => 470.0
  | _  => 0.0
  end.

(** ** 8. Electron Affinity (kJ/mol, positive = exothermic) *)

Definition electron_affinity (e : Element) : R :=
  match e with
  | H  => 72.8   | He => 0.0
  | Li => 59.6   | Be => 0.0    | B  => 26.7   | C  => 121.8
  | N  => 0.0    | O  => 141.0  | F  => 328.2  | Ne => 0.0
  | Na => 52.9   | Mg => 0.0    | Al => 41.8   | Si => 134.1
  | P  => 72.0   | S  => 200.4  | Cl => 348.6  | Ar => 0.0
  | K  => 48.4   | Ca => 2.4    | Sc => 18.1   | Ti => 7.6
  | V  => 50.7   | Cr => 64.3   | Mn => 0.0    | Fe => 15.7
  | Co => 63.7   | Ni => 112.0  | Cu => 118.4  | Zn => 0.0
  | Ga => 41.5   | Ge => 119.0  | As => 78.2   | Se => 195.0
  | Br => 324.6  | Kr => 0.0
  | Rb => 46.9   | Sr => 5.0    | Y  => 29.6   | Zr => 41.1
  | Nb => 86.1   | Mo => 71.9   | Tc => 53.0   | Ru => 101.3
  | Rh => 109.7  | Pd => 53.7   | Ag => 125.6  | Cd => 0.0
  | In => 28.9   | Sn => 107.3  | Sb => 101.1  | Te => 190.2
  | I  => 295.2  | Xe => 0.0
  | _  => 0.0
  end.

(** ** 9. Period (row in periodic table) *)

Definition period (e : Element) : nat :=
  match e with
  | H  | He => 1
  | Li | Be | B  | C  | N  | O  | F  | Ne => 2
  | Na | Mg | Al | Si | P  | S  | Cl | Ar => 3
  | K  | Ca | Sc | Ti | V  | Cr | Mn | Fe | Co | Ni | Cu | Zn
  | Ga | Ge | As | Se | Br | Kr => 4
  | Rb | Sr | Y  | Zr | Nb | Mo | Tc | Ru | Rh | Pd | Ag | Cd
  | In | Sn | Sb | Te | I  | Xe => 5
  | Cs | Ba | La | Ce | Pr | Nd | Pm | Sm | Eu | Gd | Tb | Dy
  | Ho | Er | Tm | Yb | Lu | Hf | Ta | W  | Re | Os | Ir | Pt
  | Au | Hg | Tl | Pb | Bi | Po | At | Rn => 6
  | Fr | Ra | Ac | Th | Pa | U  | Np | Pu | Am | Cm | Bk | Cf
  | Es | Fm | Md | No | Lr | Rf | Db | Sg | Bh | Hs | Mt | Ds
  | Rg | Cn | Nh | Fl | Mc | Lv | Ts | Og => 7
  end.

(** ** 10. Group (column in periodic table, using IUPAC numbering 1–18) *)

Definition group (e : Element) : nat :=
  match e with
  | H  | Li | Na | K  | Rb | Cs | Fr => 1
  | He | Ne | Ar | Kr | Xe | Rn | Og => 18
  | Be | Mg | Ca | Sr | Ba | Ra => 2
  | B  | Al | Ga | In | Tl | Nh => 13
  | C  | Si | Ge | Sn | Pb | Fl => 14
  | N  | P  | As | Sb | Bi | Mc => 15
  | O  | S  | Se | Te | Po | Lv => 16
  | F  | Cl | Br | I  | At | Ts => 17
  | Sc | Y  | La | Ac => 3
  | Ti | Zr | Hf | Rf => 4
  | V  | Nb | Ta | Db => 5
  | Cr | Mo | W  | Sg => 6
  | Mn | Tc | Re | Bh => 7
  | Fe | Ru | Os | Hs => 8
  | Co | Rh | Ir | Mt => 9
  | Ni | Pd | Pt | Ds => 10
  | Cu | Ag | Au | Rg => 11
  | Zn | Cd | Hg | Cn => 12
  (* Lanthanides and actinides: assigned to group 0 as they span 3-17 *)
  | Ce | Pr | Nd | Pm | Sm | Eu | Gd | Tb | Dy | Ho | Er | Tm | Yb | Lu => 0
  | Th | Pa | U  | Np | Pu | Am | Cm | Bk | Cf | Es | Fm | Md | No | Lr => 0
  end.

(** ** 11. Valence Electrons *)

Definition valence_electrons (e : Element) : nat :=
  match e with
  | H  => 1 | He => 2
  | Li => 1 | Be => 2 | B  => 3 | C  => 4 | N  => 5 | O  => 6 | F  => 7 | Ne => 8
  | Na => 1 | Mg => 2 | Al => 3 | Si => 4 | P  => 5 | S  => 6 | Cl => 7 | Ar => 8
  | K  => 1 | Ca => 2
  | Sc => 3 | Ti => 4 | V  => 5 | Cr => 6 | Mn => 7 | Fe => 8 | Co => 9 | Ni => 10
  | Cu => 11| Zn => 12| Ga => 3 | Ge => 4 | As => 5 | Se => 6 | Br => 7 | Kr => 8
  | Rb => 1 | Sr => 2
  | Y  => 3 | Zr => 4 | Nb => 5 | Mo => 6 | Tc => 7 | Ru => 8 | Rh => 9 | Pd => 10
  | Ag => 11| Cd => 12| In => 3 | Sn => 4 | Sb => 5 | Te => 6 | I  => 7 | Xe => 8
  | Cs => 1 | Ba => 2
  | La => 3 | Ce => 4 | Pr => 5 | Nd => 6 | Pm => 7 | Sm => 8 | Eu => 9 | Gd => 10
  | Tb => 11| Dy => 12| Ho => 13| Er => 14| Tm => 15| Yb => 16| Lu => 3
  | Hf => 4 | Ta => 5 | W  => 6 | Re => 7 | Os => 8 | Ir => 9 | Pt => 10| Au => 11
  | Hg => 12| Tl => 3 | Pb => 4 | Bi => 5 | Po => 6 | At => 7 | Rn => 8
  | Fr => 1 | Ra => 2
  | Ac => 3 | Th => 4 | Pa => 5 | U  => 6 | Np => 7 | Pu => 8 | Am => 9 | Cm => 10
  | Bk => 11| Cf => 12| Es => 13| Fm => 14| Md => 15| No => 16| Lr => 3
  | Rf => 4 | Db => 5 | Sg => 6 | Bh => 7 | Hs => 8 | Mt => 9 | Ds => 10| Rg => 11
  | Cn => 12| Nh => 3 | Fl => 4 | Mc => 5 | Lv => 6 | Ts => 7 | Og => 8
  end.

(** ** 12. Maximum Valence (typical, ignoring hypervalency) *)

Definition max_valence (e : Element) : nat :=
  match e with
  | H                     => 1
  | He | Ne | Ar | Kr     => 0
  | Li | Na | K  | Rb | Cs | Fr => 1
  | Be | Mg | Ca | Sr | Ba | Ra => 2
  | B                     => 3
  | C                     => 4
  | N                     => 3
  | O                     => 2
  | F                     => 1
  | Al                    => 3
  | Si                    => 4
  | P                     => 5   (* hypervalent: PCl₅ *)
  | S                     => 6   (* hypervalent: SF₆ *)
  | Cl                    => 7   (* ClO₄⁻ *)
  | Br                    => 5
  | I                     => 7
  | Xe                    => 8
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
  | H  | He => S_orb
  | Li | Be => S_orb
  | B  | C  | N  | O  | F  | Ne => P_orb
  | Na | Mg => S_orb
  | Al | Si | P  | S  | Cl | Ar => P_orb
  | K  | Ca => S_orb
  | Sc | Ti | V  | Cr | Mn | Fe | Co | Ni | Cu | Zn => D_orb
  | Ga | Ge | As | Se | Br | Kr => P_orb
  | Rb | Sr => S_orb
  | Y  | Zr | Nb | Mo | Tc | Ru | Rh | Pd | Ag | Cd => D_orb
  | In | Sn | Sb | Te | I  | Xe => P_orb
  | Cs | Ba => S_orb
  | La | Ce | Pr | Nd | Pm | Sm | Eu | Gd | Tb | Dy
  | Ho | Er | Tm | Yb | Lu => F_orb
  | Hf | Ta | W  | Re | Os | Ir | Pt | Au | Hg => D_orb
  | Tl | Pb | Bi | Po | At | Rn => P_orb
  | Fr | Ra => S_orb
  | Ac | Th | Pa | U  | Np | Pu | Am | Cm | Bk | Cf
  | Es | Fm | Md | No | Lr => F_orb
  | Rf | Db | Sg | Bh | Hs | Mt | Ds | Rg | Cn => D_orb
  | Nh | Fl | Mc | Lv | Ts | Og => P_orb
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
  forall e : Element, electronegativity e <= electronegativity F.
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
  intros e H.
  unfold is_hypervalent_possible in H.
  destruct (period e) eqn:Hp; try discriminate.
  destruct n; try discriminate.
  destruct n; try discriminate.
  lia.
Qed.

Close Scope R_scope.
