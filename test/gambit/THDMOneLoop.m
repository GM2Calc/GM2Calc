(*************************General 2HDM******************************)
(* Mathematica code for one-loop contributions for the general THDM 
*  to the muon anomalous magnetic moment
*
*)

(*Muon g-2 1-Loop Contributions in the Feynman Gauge*)


yffh[f_, fp_, mf_, xif_, vev_, cosab_] := 
 mf/vev*Sqrt[1 - cosab^2]*KroneckerDelta[f, fp] + 
  cosab/Sqrt[2] xif[[f, fp]]
yffH[f_, fp_, mf_, xif_, vev_, cosab_] := 
 mf/vev*cosab*KroneckerDelta[f, fp] - 
  Sqrt[1 - cosab^2]/Sqrt[2] xif[[f, fp]]
yffphi[f_, i_, j_, phi_, mf_, xif_, VCKM_, vev_, cosab_] := 
 Switch[phi, 1, yffh[i, j, mf, xif, vev, cosab], 2, 
  yffH[i, j, mf, xif, vev, cosab], 3, 
  If[f == 3, -I/Sqrt[2] xif[[i, j]], I/Sqrt[2] xif[[i, j]]], 4, 
  Switch[f, 3, -Sum[VCKM[[i, kk]]*xif[[kk, j]], {kk, 1, 3}], 2, 
   Sum[VCKM[[i, kk]]*xif[[kk, j]], {kk, 1, 3}], 1, xif[[i, j]]]]

OneLoopB[x_] := 
 If[x == 0, 2, 
  If[x == 1, 1, 
   If[x == \[Infinity], 0, (
    2 (1 - 6 x + 3 x^2 + 2 x^3 - 6 x^2 Log[x]))/(1 - x)^4]]]
OneLoopC[x_] := 
 If[x == 0, 3, 
  If[x == 1, 1, 
   If[x == \[Infinity], 0, (3 (1 - x^2 + 2 x Log[x]))/(1 - x)^3]]]
OneLoopE[x_] := 
 If[x == 0, 4, 
  If[x == 1, 1, 
   If[x == \[Infinity], 0, (
    2 (2 + 3 x - 6 x^2 + x^3 + 6 x Log[x]))/(1 - x)^4]]]
OneLoopF[x_] := 
 If[x == 0, \[Infinity], 
  If[x == 1, 1, 
   If[x == \[Infinity], 0, (3 (-3 + 4 x - x^2 - 2 Log[x]))/(
    2 (1 - x)^3)]]]

Aloop1L[f_, l_, li_, lp_, phi_, mnu_, ml_, mphi_, xiL_, VCKM_, vev_, 
  cosab_] :=
 Module[{x, term1, term2, term3},
  If[Or[phi == 1, phi == 2, phi == 3],
   x = ml[[li]]^2/mphi^2;
   term1 = ml[[l]]*ml[[l]]/ml[[l]]^2*OneLoopE[x]/24;
   term2 = ml[[l]]*ml[[lp]]/ml[[l]]^2*OneLoopE[x]/24;
   term3 = ml[[l]]*ml[[li]]/ml[[l]]^2*OneLoopF[x]/3;
   (*Print[x];Print[term1];Print[term2];Print[term3];
   Print[yffphi[f,li,lp,phi,ml[[li]],xiL,VCKM,vev,cosab]];
   Print[yffphi[f,li,l,phi,ml[[li]],xiL,VCKM,vev,cosab]];
   Print[yffphi[f,l,li,phi,ml[[l]],xiL,VCKM,vev,cosab]];
   Print[yffphi[f,lp,li,phi,ml[[lp]],xiL,VCKM,vev,cosab]];
   Print[yffphi[f,li,lp,phi,ml[[li]],xiL,VCKM,vev,cosab]];
   Print[yffphi[f,l,li,phi,ml[[l]],xiL,VCKM,vev,cosab]];*)
   Return[
    term1*
      Conjugate[
       yffphi[f, li, lp, phi, ml[[li]], xiL, VCKM, vev, cosab]]*
      yffphi[f, li, l, phi, ml[[li]], xiL, VCKM, vev, cosab] + 
     term2*Conjugate[
       yffphi[f, l, li, phi, ml[[l]], xiL, VCKM, vev, cosab]]*
      yffphi[f, lp, li, phi, ml[[lp]], xiL, VCKM, vev, cosab] + 
     term3*Conjugate[
       yffphi[f, li, lp, phi, ml[[li]], xiL, VCKM, vev, cosab]]*
      Conjugate[
       yffphi[f, l, li, phi, ml[[l]], xiL, VCKM, vev, cosab]]],
   If[phi == 4,
    x = mnu[[li]]^2/mphi^2;
    term1 = ml[[l]]/ml[[l]]*OneLoopB[x]/24;
    (*Print[x];Print[term1];
    Print[yffphi[f,li,lp,phi,ml[[li]],xiL,VCKM,vev,cosab]];
    Print[yffphi[f,li,l,phi,ml[[li]],xiL,VCKM,vev,cosab]];*)
    term1*
     Conjugate[
      yffphi[f, li, lp, phi, ml[[li]], xiL, VCKM, vev, cosab]]*
     yffphi[f, li, l, phi, ml[[li]], xiL, VCKM, vev, cosab]
    ]
   ]
  ]
Aloop1R[f_, l_, li_, lp_, phi_, mnu_, ml_, mphi_, xiL_, VCKM_, vev_, 
  cosab_] :=
 Module[{x, term1, term2, term3},
  If[Or[phi == 1, phi == 2, phi == 3],
   x = ml[[li]]^2/mphi^2;
   term1 = ml[[l]]*ml[[l]]/ml[[l]]^2*OneLoopE[x]/24;
   term2 = ml[[l]]*ml[[lp]]/ml[[l]]^2*OneLoopE[x]/24;
   term3 = ml[[l]]*ml[[li]]/ml[[l]]^2*OneLoopF[x]/3;
   term1*Conjugate[
      yffphi[f, l, li, phi, ml[[l]], xiL, VCKM, vev, cosab]]*
     yffphi[f, lp, li, phi, ml[[lp]], xiL, VCKM, vev, cosab] + 
    term2*Conjugate[
      yffphi[f, li, lp, phi, ml[[li]], xiL, VCKM, vev, cosab]]*
     yffphi[f, li, l, phi, ml[[li]], xiL, VCKM, vev, cosab] + 
    term3*Conjugate[
      yffphi[f, li, l, phi, ml[[li]], xiL, VCKM, vev, cosab]]*
     Conjugate[
      yffphi[f, lp, li, phi, ml[[lp]], xiL, VCKM, vev, cosab]],
   If[phi == 4,
    x = mnu[[l]]^2/mphi^2;
    term1 = ml[[lp]]/ml[[l]]*OneLoopB[x]/24;
    term1*
     Conjugate[
      yffphi[f, li, lp, phi, ml[[li]], xiL, VCKM, vev, cosab]]*
     yffphi[f, li, l, phi, ml[[li]], xiL, VCKM, vev, cosab]
    ]
   ]
  ]
