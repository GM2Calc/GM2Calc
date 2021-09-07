(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

Off[General::spell]

Model`Name = "GTHDM";
Model`NameLaTeX ="General Two Higgs Doublet Model";
Model`Authors = "D.Jacob,F.Staub,A.Voigt";
Model`Date = "2021-07-07";

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Gauge Superfields *)

Gauge[[1]]={B,   U[1], hypercharge, g1, False};
Gauge[[2]]={WB, SU[2], left,        g2, True };
Gauge[[3]]={G,  SU[3], color,       g3, False};


(* Chiral Superfields *)

FermionFields[[1]] = {q, 3, {uL,dL},     1/6, 2,  3};
FermionFields[[2]] = {l, 3, {vL,eL},    -1/2, 2,  1};
FermionFields[[3]] = {d, 3, conj[dR],    1/3, 1, -3};
FermionFields[[4]] = {u, 3, conj[uR],   -2/3, 1, -3};
FermionFields[[5]] = {e, 3, conj[eR],      1, 1,  1};

ScalarFields[[1]]  = {H1, 1, {H1p, H10}, 1/2, 2,  1};
ScalarFields[[2]]  = {H2, 1, {H2p, H20}, 1/2, 2,  1};

(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][Additional]= {
	{LagHC, { AddHC->True}},
	{LagNoHC,{ AddHC->False}}
};

LagNoHC = -(+ M112 conj[H1].H1
            + M222 conj[H2].H2
            + Lambda1/2 conj[H1].H1.conj[H1].H1
            + Lambda2/2 conj[H2].H2.conj[H2].H2
            + Lambda3 conj[H2].H2.conj[H1].H1
            + Lambda4 conj[H2].H1.conj[H1].H2);

LagHC = -(-M122 conj[H1].H2
          + Lambda5/2 conj[H2].H1.conj[H2].H1
          + Lambda6 conj[H1].H1.conj[H1].H2
          + Lambda7 conj[H2].H2.conj[H1].H2
          + Gammad q.conj[H1].d + Gammae l.conj[H1].e + Piu q.H2.u
          + Pid q.conj[H2].d + Pie l.conj[H2].e + Gammau q.H1.u);

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] = {
    {{VB,VWB[3]},{VP,VZ},ZZ},
    {{VWB[1],VWB[2]},{VWm,conj[VWm]},ZW}
};

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs] = {
    {H10, {v1, 1/Sqrt[2]}, {sigma1, \[ImaginaryI]/Sqrt[2]},{phi1, 1/Sqrt[2]}},
    {H20, {v2, 1/Sqrt[2]}, {sigma2, \[ImaginaryI]/Sqrt[2]},{phi2, 1/Sqrt[2]}}
};


DEFINITION[EWSB][MatterSector] = {
    {{phi1, phi2}, {hh, ZH}},
    {{sigma1, sigma2}, {Ah, ZA}},
    {{conj[H1p],conj[H2p]},{Hm,ZP}},
    {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
    {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
    {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}}
};

(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][DiracSpinors] = {
    Fd ->{DL, conj[DR]},
    Fe ->{EL, conj[ER]},
    Fu ->{UL, conj[UR]},
    Fv ->{vL, 0}
};

DEFINITION[EWSB][GaugeES] = {
    Fd1 ->{FdL, 0},
    Fd2 ->{0, FdR},
    Fu1 ->{Fu1, 0},
    Fu2 ->{0, Fu2},
    Fe1 ->{Fe1, 0},
    Fe2 ->{0, Fe2}
};
