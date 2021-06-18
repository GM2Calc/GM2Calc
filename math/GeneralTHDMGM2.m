(* ::Package:: *)

(* ::Title:: *)
(*General 2HDM*)


(* ::Chapter:: *)
(*Muon g-2 2-Loop in the Feynman Gauge*)


(* ::Section:: *)
(*Barr-Zee (2 boson) 2-loop Contributions*)


(* ::Text:: *)
(*(*All functions taken from arXiv:1502.04199v3*)*)


(* ::Text:: *)
(*(*Equations (25-29)*)*)


(* ::Input:: *)
(*F1[w_]:=Module[{x},w/2 NIntegrate[(2x(1-x)-1)/(w-x(1-x)) Log[w/(x(1-x))],{x,0,1}]]*)
(*FF1[w_]:=Module[{x},w/2 NIntegrate[1/(w-x(1-x)) Log[w/(x(1-x))],{x,0,1}]]*)
(*F2[w_]:=Module[{x},1/2 NIntegrate[(x(x-1))/(w-x(1-x)) Log[w/(x(1-x))],{x,0,1}]]*)
(*F3[w_]:=Module[{x},1/2 NIntegrate[(x w(3x(4x-1)+10)-x(1-x))/(w-x(1-x)) Log[w/(x(1-x))],{x,0,1}]]*)
(*G[wa_,wb_,x_]:=Log[(wa*x+wb*(1-x))/(x(1-x))]/(x(1-x)-wa*x-wb*(1-x))*)


(* ::Text:: *)
(*(*Equations (19-21)*)*)


(* ::Input:: *)
(*BarrZee\[Phi]\[Gamma]F[MS_,mf_,\[Alpha]_,m\[Mu]_,v_,Qf_,Nc_,YfS_,Y\[Mu]S_]:=(\[Alpha] m\[Mu]^2)/(4\[Pi]^3v^2) Nc (Qf)^2(Re[YfS]Re[Y\[Mu]S]F1[mf^2/MS^2]+Im[YfS]Im[Y\[Mu]S]FF1[mf^2/MS^2])*)
(*BarrZee\[Phi]\[Gamma]C[MS_,MC_,\[Alpha]_,m\[Mu]_,v_,\[Lambda]\[Phi]CC_,YlS_]:=(\[Alpha] m\[Mu]^2)/(8\[Pi]^3MS^2) Re[YlS]\[Lambda]\[Phi]CC F2[MC^2/MS^2]*)
(*BarrZee\[Phi]\[Gamma]W[MS_,\[Alpha]_,m\[Mu]_,MW_,v_,g\[Phi]WW_,YlS_]:=(\[Alpha] m\[Mu]^2)/(8\[Pi]^3v^2) Re[YlS]g\[Phi]WW F3[MW^2/MS^2]*)


(* ::Text:: *)
(*(*Equations (22-24)*)*)


(* ::Input:: *)
(*BarrZeeCWF[MC_,mf_,\[Alpha]_,m\[Mu]_,MW_,sw2_,v_,Qf_,Nc_,VCKM_,YfA_,Y\[Mu]A_]:=(\[Alpha] m\[Mu]^2Nc Abs[VCKM[[3,2]]]^2)/(32\[Pi]^3sw2 v^2(MC^2-MW^2)) NIntegrate[(Qf[[3]] x+Qf[[2]](1-x))*(Re[YfA[[2]] Conjugate[Y\[Mu]A]]mf[[2]]^2 x(1-x)+Re[YfA[[3]] Conjugate[Y\[Mu]A]]mf[[3]]^2x(1+x))*(G[mf[[3]]^2/MC^2,mf[[2]]^2/MC^2,x]-G[mf[[3]]^2/MW^2,mf[[2]]^2/MW^2,x]),{x,0,1}]*)
(*BarrZeeCWC[MS_,MC_,\[Alpha]_,m\[Mu]_,MW_,sw2_,v_,\[Lambda]\[Phi]CC_,\[Lambda]\[Phi]CW_,YlA_]:=(\[Alpha] m\[Mu]^2)/(64\[Pi]^3sw2 (MC^2-MW^2)) Re[Conjugate[YlA]\[Lambda]\[Phi]CW]\[Lambda]\[Phi]CC NIntegrate[x^2(x-1)*(G[1,MS^2/MC^2,x]-G[MC^2/MW^2,MS^2/MW^2,x]),{x,0,1}]*)
(*BarrZeeCWW[MS_,MC_,\[Alpha]_,m\[Mu]_,MW_,sw2_,v_,g\[Phi]WW_,\[Lambda]\[Phi]CW_,YlA_]:=(\[Alpha] m\[Mu]^2)/(64\[Pi]^3sw2 v^2(MC^2-MW^2)) Re[Conjugate[YlA]g\[Phi]WW \[Lambda]\[Phi]CW]NIntegrate[x^2*((MC^2+MW^2-MS^2)(1-x)-4MW^2)*(G[MW^2/MC^2,MS^2/MC^2,x]-G[1,MS^2/MW^2,x]),{x,0,1}]*)


(* ::Section:: *)
(*Bosonic (3 boson) 2-loop Contributions*)


(* ::Text:: *)
(*(*All functions taken from arXiv:1607.06292v2 .  Not Implemented*)*)


(* ::Text:: *)
(*(*Equation (43)*)*)


(* ::Input:: *)
(*B3WBoson[MS_,\[Alpha]_,m\[Mu]_,MW_,MZ_,sw2_,\[Epsilon]_,Cs_,YlS_]:=Module[{yS},*)
(*yS=MS^2/MW^2;*)
(*Return[(Cs YlS \[Alpha]^2)/(576\[Pi]^2(1-sw2)^2sw2^2) m\[Mu]^2/MZ^2 (3/\[Epsilon]-6L[MW^2]-55/2+32/yS-(4\[Pi]^2)/3 (4+3yS)/yS^2-(35+32/yS)Log[yS]+(6+32/yS^2+24/yS-32yS)PolyLog[1-yS,2]+((10+70yS-32yS^2)/((yS-4)yS))\[CapitalPhi][Sqrt[yS],1,1])]]*)


(* ::Text:: *)
(*(*Equations (45-46)*)*)


(* ::Input:: *)
(*fa[x_]:=(3(4-x))/x-\[Pi]^2/2 (4+3x)/x^2-(3(4+x))/x Log[x]+(12+9x-3x^3)/x^2 PolyLog[1-x,2]+(3(2+x))/x \[CapitalPhi][Sqrt[x],1,1]*)
(*fb[x_]:=(\[Pi]^2(8+6x-12x^4+3x^5))/x^2+(6(-8+2x+3x^2))/x+(12(4+x+3x^2))/x Log[x]+9(-4+x)x^2Log[x]^2+(12(-4-3x+4x^3-12x^4+3x^5))/x^2 PolyLog[1-x,2]+(6(4+2x-6x^2+3x^3))/x \[CapitalPhi][Sqrt[x],1,1]*)


(* ::Text:: *)
(*(*Equation (44)*)*)


(* ::Input:: *)
(*B3ZBoson[MS_,\[Alpha]_,m\[Mu]_,MW_,MZ_,sw2_,Cs_,YlS_]:=(Cs YlS \[Alpha]^2)/(576\[Pi]^2(1-sw2)^2sw2^2) m\[Mu]^2/MZ^2 (fa[MS^2/MZ^2]+sw2(1-2sw2)fb[MS^2/MZ^2])*)


(* ::Section:: *)
(*Fermionic 2-loop Contributions*)


(* ::Text:: *)
(*(*All functions taken from arXiv:1607.06292v2*)*)


(* ::Text:: *)
(*(*Equations (68-70)*)*)


(* ::Input:: *)
(*\[CapitalPhi][mm1_,mm2_,mm3_]:=Module[{\[Lambda],\[Alpha]p,\[Alpha]m,m1,m2,m3},*)
(*{m1,m2,m3}=Sort[{mm1,mm2,mm3}];*)
(*\[Lambda]=Sqrt[m1^4+m2^4+m3^4-2m1^2m2^2-2m2^2m3^2-2m3^2m1^2];*)
(*\[Alpha]p=(m3^2+m1^2-m2^2-\[Lambda])/(2m3^2);*)
(*\[Alpha]m=(m3^2-m1^2+m2^2-\[Lambda])/(2m3^2);*)
(*Return[\[Lambda]/2 (2Log[\[Alpha]p]Log[\[Alpha]m]-Log[m1^2/m3^2]Log[m2^2/m3^2]-2PolyLog[2,\[Alpha]p]-2PolyLog[2,\[Alpha]m]+\[Pi]^2/3)]*)
(*]*)


(* ::Text:: *)
(*(*Equations (56-57)*)*)


(* ::Input:: *)
(*F\[Phi][MS_,mf_]:=-2+Log[MS^2/mf^2]-((MS^2-2mf^2)/MS^2) \[CapitalPhi][MS,mf,mf]/(MS^2-4mf^2)*)
(*FA[MS_,mf_]:=\[CapitalPhi][MS,mf,mf]/(MS^2-4mf^2)*)


(* ::Text:: *)
(*(*Equations (54-55)*)*)


(* ::Input:: *)
(*f\[Gamma]\[Phi][MS_,mf_,\[Alpha]_,m\[Mu]_,MW_,sw2_,Qf_,Nc_]:=(\[Alpha]^2m\[Mu]^2)/(4\[Pi]^2MW^2sw2) (Qf^2Nc)(mf^2/MS^2)F\[Phi][MS,mf]*)
(*fZ\[Phi][MS_,mf_,\[Alpha]_,m\[Mu]_,MW_,MZ_,sw2_,Qf_,Nc_,glv_,gfv_]:=(\[Alpha]^2m\[Mu]^2)/(4\[Pi]^2MW^2sw2) (- ((Qf Nc glv gfv)/(sw2(1-sw2))))(mf^2/(MS^2-MZ^2))(F\[Phi][MS,mf]-F\[Phi][MZ,mf])*)
(*f\[Gamma]A[MA_,mf_,\[Alpha]_,m\[Mu]_,MW_,sw2_,Qf_,Nc_]:=(\[Alpha]^2m\[Mu]^2)/(4\[Pi]^2MW^2sw2) (Qf^2Nc)(mf^2/MA^2)FA[MA,mf]*)
(*fZA[MA_,mf_,\[Alpha]_,m\[Mu]_,MW_,MZ_,sw2_,Qf_,Nc_,glv_,gfv_]:=(\[Alpha]^2m\[Mu]^2)/(4\[Pi]^2MW^2sw2) (- ((Qf Nc glv gfv)/(sw2(1-sw2))))(mf^2/(MA^2-MZ^2))(FA[MA,mf]-FA[MZ,mf])*)


(* ::Text:: *)
(*(*Equations (53)*)*)


(* ::Input:: *)
(*FNeutral[Mh_,MH_,MA_,ml_,md_,mu_,\[Alpha]_,m\[Mu]_,MW_,MZ_,sw2_,Ql_,Qd_,Qu_,Nc_,glv_,gdv_,guv_,Ylh_,Ydh_,Yuh_,YlH_,YdH_,YuH_,YlA_,YdA_,YuA_]:=(f\[Gamma]\[Phi][Mh,ml,\[Alpha],m\[Mu],MW,sw2,Ql,1]+fZ\[Phi][Mh,ml,\[Alpha],m\[Mu],MW,MZ,sw2,Ql,1,glv,glv])Ylh Ylh+(f\[Gamma]\[Phi][Mh,md,\[Alpha],m\[Mu],MW,sw2,Qd,Nc]+fZ\[Phi][Mh,md,\[Alpha],m\[Mu],MW,MZ,sw2,Qd,Nc,glv,gdv])Ydh Ylh+(f\[Gamma]\[Phi][Mh,mu,\[Alpha],m\[Mu],MW,sw2,Qu,Nc]+fZ\[Phi][Mh,mu,\[Alpha],m\[Mu],MW,MZ,sw2,Qu,Nc,glv,guv])Yuh Ylh+(f\[Gamma]\[Phi][MH,ml,\[Alpha],m\[Mu],MW,sw2,Ql,1]+fZ\[Phi][MH,ml,\[Alpha],m\[Mu],MW,MZ,sw2,Ql,1,glv,glv])YlH YlH+(f\[Gamma]\[Phi][MH,md,\[Alpha],m\[Mu],MW,sw2,Qd,Nc]+fZ\[Phi][MH,md,\[Alpha],m\[Mu],MW,MZ,sw2,Qd,Nc,glv,gdv])YdH YlH+(f\[Gamma]\[Phi][MH,mu,\[Alpha],m\[Mu],MW,sw2,Qu,Nc]+fZ\[Phi][MH,mu,\[Alpha],m\[Mu],MW,MZ,sw2,Qu,Nc,glv,guv])YuH YlH+(f\[Gamma]\[Phi][MA,ml,\[Alpha],m\[Mu],MW,sw2,Ql,1]+fZ\[Phi][MA,ml,\[Alpha],m\[Mu],MW,MZ,sw2,Ql,1,glv,glv])YlA YlA+(f\[Gamma]\[Phi][MA,md,\[Alpha],m\[Mu],MW,sw2,Qd,Nc]+fZ\[Phi][MA,md,\[Alpha],m\[Mu],MW,MZ,sw2,Qd,Nc,glv,gdv])YdA YlA+(f\[Gamma]\[Phi][MA,mu,\[Alpha],m\[Mu],MW,sw2,Qu,Nc]+fZ\[Phi][MA,mu,\[Alpha],m\[Mu],MW,MZ,sw2,Qu,Nc,glv,guv])YuA YlA*)


(* ::Input:: *)
(*FNeutral[M\[Phi]_,mf_,\[Alpha]_,m\[Mu]_,MW_,MZ_,sw2_,Qf_,gfv_,Yf\[Phi]_,Y\[Mu]\[Phi]_]:=*)
(*Sum[*)
(*(f\[Gamma]\[Phi][M\[Phi][[i]],mf[[j]],\[Alpha],m\[Mu],MW,sw2,Qf[[j]],{1,3,3}[[j]]]+fZ\[Phi][M\[Phi][[i]],mf[[j]],\[Alpha],m\[Mu],MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]])Yf\[Phi][[i]][[j]] Y\[Mu]\[Phi][[i]],*)
(*{j,1,3},*)
(*{i,1,3}]*)


(* ::Text:: *)
(*(*Equations (60-62)*)*)


(* ::Input:: *)
(*FCl[xl_]:=xl+xl(xl-1)(PolyLog[2,1-1/xl]-\[Pi]^2/6)+(xl-1/2)Log[xl]*)
(*FCd[xu_,xd_,Qu_,Qd_]:=Module[{c,cb,y,s},*)
(*c=(xu-xd)^2-Qu xu+Qd xd;*)
(*cb=(xu-Qu)xu-(xd+Qd)xd;*)
(*y=(xu-xd)^2-2(xu+xd)+1;*)
(*s=(Qu+Qd)/4;*)
(*Return[-(xu-xd)+(cb/y-c((xu-xd)/y))\[CapitalPhi][Sqrt[xd],Sqrt[xu],1]+c(PolyLog[2,1-xd/xu]-1/2 Log[xu]Log[xd/xu](*\[CapitalPhi][Sqrt[xd],Sqrt[xu],1]*))+(s+xd)Log[xd]+(s-xu)Log[xu]]*)
(*]*)
(*FCu[xu_,xd_,Qu_,Qd_]:=Module[{c,cb,y,s},*)
(*c=(xu-xd)^2-Qu xu+Qd xd;*)
(*cb=(xu-Qu)xu-(xd+Qd)xd;*)
(*y=(xu-xd)^2-2(xu+xd)+1;*)
(*s=(Qu+Qd)/4;*)
(*Return[-(xu-xd)+(cb/y-c((xu-xd)/y))\[CapitalPhi][Sqrt[xd],Sqrt[xu],1]+c(PolyLog[2,1-xd/xu]-1/2 Log[xu]Log[xd/xu](*\[CapitalPhi][Sqrt[xd],Sqrt[xu],1]*))+(s+xd)Log[xd]+(s-xu)Log[xu]-4/3 ((xu-xd-1)/y)\[CapitalPhi][Sqrt[xd],Sqrt[xu],1]-1/3 (Log[xd]^2-Log[xu]^2)]]*)


(* ::Text:: *)
(*(*Equation (59)*)*)


(* ::Input:: *)
(*fCl[MC_,ml_,\[Alpha]_,m\[Mu]_,MW_,sw2_]:=(\[Alpha]^2m\[Mu]^2)/(32\[Pi]^2MW^2sw2^2) ml^2/(MC^2-MW^2) (FCl[ml^2/MC^2]-FCl[ml^2/MW^2])*)
(*fCd[MC_,mu_,md_,\[Alpha]_,m\[Mu]_,MW_,sw2_,Qu_,Qd_,Nc_]:=(\[Alpha]^2m\[Mu]^2)/(32\[Pi]^2MW^2sw2^2) (Nc md^2)/(MC^2-MW^2) (FCd[mu^2/MC^2,md^2/MC^2,Qu,Qd]-FCd[mu^2/MW^2,md^2/MW^2,Qu,Qd])*)
(*fCu[MC_,mu_,md_,\[Alpha]_,m\[Mu]_,MW_,sw2_,Qu_,Qd_,Nc_]:=(\[Alpha]^2m\[Mu]^2)/(32\[Pi]^2MW^2sw2^2) (Nc mu^2)/(MC^2-MW^2) (FCu[mu^2/MC^2,md^2/MC^2,Qu+2,Qd+2]-FCu[mu^2/MW^2,md^2/MW^2,Qu+2,Qd+2])*)


(* ::Text:: *)
(*(*Equation (58)*)*)


(* ::Input:: *)
(*FCharged[MC_,mf_,\[Alpha]_,m\[Mu]_,MW_,sw2_,Qf_,YfA_,Y\[Mu]A_]:=fCl[MC,mf[[1]],\[Alpha],m\[Mu],MW,sw2]YfA[[1]] Y\[Mu]A+fCd[MC,mf[[3]],mf[[2]],\[Alpha],m\[Mu],MW,sw2,Qf[[3]],Qf[[2]],3]YfA[[2]] Y\[Mu]A+fCu[MC,mf[[3]],mf[[2]],\[Alpha],m\[Mu],MW,sw2,Qf[[3]],Qf[[2]],3]YfA[[3]] Y\[Mu]A*)


(* ::Input:: *)
(*Fermionic2Loop[M\[CapitalPhi]_,MhSM_,mf_,\[Alpha]_,m\[Mu]_,MW_,MZ_,sw2_,Qf_,gfv_,YfS_,Y\[Mu]S_,YfhSM_]:=FNeutral[M\[CapitalPhi][[1]],M\[CapitalPhi][[2]],M\[CapitalPhi][[3]],mf,\[Alpha],m\[Mu],MW,MZ,sw2,Qf,gfv,YfS,Y\[Mu]S]+*)
(*FCharged[M\[CapitalPhi][[4]],mf,\[Alpha],m\[Mu],MW,sw2,Qf,I*YfS[[3,;;]],I*Y\[Mu]S[[3]]]-*)
(*Sum[( *)
(*f\[Gamma]\[Phi][MhSM,mf[[j]],\[Alpha],m\[Mu],MW,sw2,Qf[[j]],{1,3,3}[[j]]]+fZ\[Phi][MhSM,mf[[j]],\[Alpha],m\[Mu],MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]]*)
(*)(*YfhSM\[LeftDoubleBracket]j\[RightDoubleBracket] YfhSM\[LeftDoubleBracket]1\[RightDoubleBracket]*),{j,1,3}]*)
