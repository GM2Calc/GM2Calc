(*************************General 2HDM******************************)
(* Mathematica code for two-loop contributions for the general THDM 
*  to the muon anomalous magnetic moment
*
*  Comments: [1] The 3-boson 2 loop contributions have not been 
*                 verified with arXiv:1502.04199v3, nor been 
*                 implemented into a .cpp file
*)

(*Muon g-2 2-Loop Contributions in the Feynman Gauge*)


(*----Barr-Zee (2 boson) 2-loop Contributions----*)


(*(*All functions taken from arXiv:1502.04199v3*)*)


(*(*Equations (25-29)*)*)


F1[w_]:=Module[{x},w/2 NIntegrate[(2x(1-x)-1)/(w-x(1-x)) Log[w/(x(1-x))],{x,0,1}]]
FF1[w_]:=Module[{x},w/2 NIntegrate[1/(w-x(1-x)) Log[w/(x(1-x))],{x,0,1}]]
F2[w_]:=Module[{x},1/2 NIntegrate[(x(x-1))/(w-x(1-x)) Log[w/(x(1-x))],{x,0,1}]]
F3[w_]:=Module[{x},1/2 NIntegrate[(x w(3x(4x-1)+10)-x(1-x))/(w-x(1-x)) Log[w/(x(1-x))],{x,0,1}]]
G[wa_,wb_,x_]:=Log[(wa*x+wb*(1-x))/(x(1-x))]/(x(1-x)-wa*x-wb*(1-x))


(*(*Equations (19-21)*)*)


BarrZeephigammaF[MS_,mf_,alpha_,mmu_,v_,Qf_,Nc_,YfS_,YmuS_]:=(alpha mmu^2)/(4Pi^3v^2) Nc (Qf)^2(Re[YfS]Re[YmuS]F1[mf^2/MS^2]+Im[YfS]Im[YmuS]FF1[mf^2/MS^2])
BarrZeephigammaC[MS_,MC_,alpha_,mmu_,v_,lambdaphiCC_,YlS_]:=(alpha mmu^2)/(8Pi^3MS^2) Re[YlS]lambdaphiCC F2[MC^2/MS^2]
BarrZeephigammaW[MS_,alpha_,mmu_,MW_,v_,gphiWW_,YlS_]:=(alpha mmu^2)/(8Pi^3v^2) Re[YlS]gphiWW F3[MW^2/MS^2]


(*(*Equations (22-24)*)*)


BarrZeeCWF[MC_,mf_,alpha_,mmu_,MW_,sw2_,v_,Qf_,Nc_,VCKM_,YfA_,YmuA_]:=(alpha mmu^2Nc Abs[VCKM[[3,2]]]^2)/(32Pi^3sw2 v^2(MC^2-MW^2)) NIntegrate[(Qf[[3]] x+Qf[[2]](1-x))*(Re[YfA[[2]] Conjugate[YmuA]]mf[[2]]^2 x(1-x)+Re[YfA[[3]] Conjugate[YmuA]]mf[[3]]^2x(1+x))*(G[mf[[3]]^2/MC^2,mf[[2]]^2/MC^2,x]-G[mf[[3]]^2/MW^2,mf[[2]]^2/MW^2,x]),{x,0,1}]
BarrZeeCWC[MS_,MC_,alpha_,mmu_,MW_,sw2_,v_,lambdaphiCC_,lambdaphiCW_,YlA_]:=(alpha mmu^2)/(64Pi^3sw2 (MC^2-MW^2)) Re[Conjugate[YlA]lambdaphiCW]lambdaphiCC NIntegrate[x^2(x-1)*(G[1,MS^2/MC^2,x]-G[MC^2/MW^2,MS^2/MW^2,x]),{x,0,1}]
BarrZeeCWW[MS_,MC_,alpha_,mmu_,MW_,sw2_,v_,gphiWW_,lambdaphiCW_,YlA_]:=(alpha mmu^2)/(64Pi^3sw2 v^2(MC^2-MW^2)) Re[Conjugate[YlA]gphiWW lambdaphiCW]NIntegrate[x^2*((MC^2+MW^2-MS^2)(1-x)-4MW^2)*(G[MW^2/MC^2,MS^2/MC^2,x]-G[1,MS^2/MW^2,x]),{x,0,1}]


(*----Bosonic (3 boson) 2-loop Contributions----*)


(*(*All functions taken from arXiv:1607.06292v2 .  Not Implemented in cpp file*)*)


(*(*Equation (43)*)*)


B3WBoson[MS_,alpha_,mmu_,MW_,MZ_,sw2_,\[Epsilon]_,Cs_,YlS_]:=Module[{yS},
yS=MS^2/MW^2;
Return[(Cs YlS alpha^2)/(576Pi^2(1-sw2)^2sw2^2) mmu^2/MZ^2 (3/\[Epsilon]-6L[MW^2]-55/2+32/yS-(4Pi^2)/3 (4+3yS)/yS^2-(35+32/yS)Log[yS]+(6+32/yS^2+24/yS-32yS)PolyLog[1-yS,2]+((10+70yS-32yS^2)/((yS-4)yS))Phi[Sqrt[yS],1,1])]]


(*(*Equations (45-46)*)*)


fa[x_]:=(3(4-x))/x-Pi^2/2 (4+3x)/x^2-(3(4+x))/x Log[x]+(12+9x-3x^3)/x^2 PolyLog[1-x,2]+(3(2+x))/x Phi[Sqrt[x],1,1]
fb[x_]:=(Pi^2(8+6x-12x^4+3x^5))/x^2+(6(-8+2x+3x^2))/x+(12(4+x+3x^2))/x Log[x]+9(-4+x)x^2Log[x]^2+(12(-4-3x+4x^3-12x^4+3x^5))/x^2 PolyLog[1-x,2]+(6(4+2x-6x^2+3x^3))/x Phi[Sqrt[x],1,1]


(*(*Equation (44)*)*)


B3ZBoson[MS_,alpha_,mmu_,MW_,MZ_,sw2_,Cs_,YlS_]:=(Cs YlS alpha^2)/(576Pi^2(1-sw2)^2sw2^2) mmu^2/MZ^2 (fa[MS^2/MZ^2]+sw2(1-2sw2)fb[MS^2/MZ^2])


(*--------------------------------------*)
(*----Fermionic 2-loop Contributions----*)
(*--------------------------------------*)

(*All functions taken from arXiv:1607.06292v2*)


(*Equations (68-70)*)


Phi[mm1_,mm2_,mm3_]:=Module[{lambda,alphap,alpham,m1,m2,m3},
{m1,m2,m3}=Sort[{mm1,mm2,mm3}];
lambda=Sqrt[m1^4+m2^4+m3^4-2m1^2m2^2-2m2^2m3^2-2m3^2m1^2];
alphap=(m3^2+m1^2-m2^2-lambda)/(2m3^2);
alpham=(m3^2-m1^2+m2^2-lambda)/(2m3^2);
Return[lambda/2 (2Log[alphap]Log[alpham]-Log[m1^2/m3^2]Log[m2^2/m3^2]-2PolyLog[2,alphap]-2PolyLog[2,alpham]+Pi^2/3)]
]


(*Equations (56-57)*)


Fphi[MS_,mf_]:=-2+Log[MS^2/mf^2]-((MS^2-2mf^2)/MS^2) Phi[MS,mf,mf]/(MS^2-4mf^2)
FA[MS_,mf_]:=Phi[MS,mf,mf]/(MS^2-4mf^2)


(*Equations (54-55)*)


fgammaphi[MS_,mf_,alpha_,mmu_,MW_,sw2_,Qf_,Nc_]:=(alpha^2mmu^2)/(4Pi^2MW^2sw2) (Qf^2Nc)(mf^2/MS^2)Fphi[MS,mf]
fZphi[MS_,mf_,alpha_,mmu_,MW_,MZ_,sw2_,Qf_,Nc_,glv_,gfv_]:=(alpha^2mmu^2)/(4Pi^2MW^2sw2) (- ((Qf Nc glv gfv)/(sw2(1-sw2))))(mf^2/(MS^2-MZ^2))(Fphi[MS,mf]-Fphi[MZ,mf])
fgammaA[MA_,mf_,alpha_,mmu_,MW_,sw2_,Qf_,Nc_]:=(alpha^2mmu^2)/(4Pi^2MW^2sw2) (Qf^2Nc)(mf^2/MA^2)FA[MA,mf]
fZA[MA_,mf_,alpha_,mmu_,MW_,MZ_,sw2_,Qf_,Nc_,glv_,gfv_]:=(alpha^2mmu^2)/(4Pi^2MW^2sw2) (- ((Qf Nc glv gfv)/(sw2(1-sw2))))(mf^2/(MA^2-MZ^2))(FA[MA,mf]-FA[MZ,mf])


(*Equations (53)*)


FNeutral[Mh_,MH_,MA_,ml_,md_,mu_,alpha_,mmu_,MW_,MZ_,sw2_,Ql_,Qd_,Qu_,Nc_,glv_,gdv_,guv_,Ylh_,Ydh_,Yuh_,YlH_,YdH_,YuH_,YlA_,YdA_,YuA_]:=(fgammaphi[Mh,ml,alpha,mmu,MW,sw2,Ql,1]+fZphi[Mh,ml,alpha,mmu,MW,MZ,sw2,Ql,1,glv,glv])Ylh Ylh+(fgammaphi[Mh,md,alpha,mmu,MW,sw2,Qd,Nc]+fZphi[Mh,md,alpha,mmu,MW,MZ,sw2,Qd,Nc,glv,gdv])Ydh Ylh+(fgammaphi[Mh,mu,alpha,mmu,MW,sw2,Qu,Nc]+fZphi[Mh,mu,alpha,mmu,MW,MZ,sw2,Qu,Nc,glv,guv])Yuh Ylh+(fgammaphi[MH,ml,alpha,mmu,MW,sw2,Ql,1]+fZphi[MH,ml,alpha,mmu,MW,MZ,sw2,Ql,1,glv,glv])YlH YlH+(fgammaphi[MH,md,alpha,mmu,MW,sw2,Qd,Nc]+fZphi[MH,md,alpha,mmu,MW,MZ,sw2,Qd,Nc,glv,gdv])YdH YlH+(fgammaphi[MH,mu,alpha,mmu,MW,sw2,Qu,Nc]+fZphi[MH,mu,alpha,mmu,MW,MZ,sw2,Qu,Nc,glv,guv])YuH YlH+(fgammaphi[MA,ml,alpha,mmu,MW,sw2,Ql,1]+fZphi[MA,ml,alpha,mmu,MW,MZ,sw2,Ql,1,glv,glv])YlA YlA+(fgammaphi[MA,md,alpha,mmu,MW,sw2,Qd,Nc]+fZphi[MA,md,alpha,mmu,MW,MZ,sw2,Qd,Nc,glv,gdv])YdA YlA+(fgammaphi[MA,mu,alpha,mmu,MW,sw2,Qu,Nc]+fZphi[MA,mu,alpha,mmu,MW,MZ,sw2,Qu,Nc,glv,guv])YuA YlA


FNeutral[Mphi_,mf_,alpha_,mmu_,MW_,MZ_,sw2_,Qf_,gfv_,Yfphi_,Ymuphi_]:=
Sum[
(fgammaphi[Mphi[[i]],mf[[j]],alpha,mmu,MW,sw2,Qf[[j]],{1,3,3}[[j]]]+fZphi[Mphi[[i]],mf[[j]],alpha,mmu,MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]])Yfphi[[i]][[j]] Ymuphi[[i]],
{j,1,3},
{i,1,3}]


(*Equations (60-62)*)


FCl[xl_]:=xl+xl(xl-1)(PolyLog[2,1-1/xl]-Pi^2/6)+(xl-1/2)Log[xl]
FCd[xu_,xd_,Qu_,Qd_]:=Module[{c,cb,y,s},
c=(xu-xd)^2-Qu xu+Qd xd;
cb=(xu-Qu)xu-(xd+Qd)xd;
y=(xu-xd)^2-2(xu+xd)+1;
s=(Qu+Qd)/4;
Return[-(xu-xd)+(cb/y-c((xu-xd)/y))Phi[Sqrt[xd],Sqrt[xu],1]+c(PolyLog[2,1-xd/xu]-1/2 Log[xu]Log[xd/xu](*Phi[Sqrt[xd],Sqrt[xu],1]*))+(s+xd)Log[xd]+(s-xu)Log[xu]]
]
FCu[xu_,xd_,Qu_,Qd_]:=Module[{c,cb,y,s},
c=(xu-xd)^2-Qu xu+Qd xd;
cb=(xu-Qu)xu-(xd+Qd)xd;
y=(xu-xd)^2-2(xu+xd)+1;
s=(Qu+Qd)/4;
Return[-(xu-xd)+(cb/y-c((xu-xd)/y))Phi[Sqrt[xd],Sqrt[xu],1]+c(PolyLog[2,1-xd/xu]-1/2 Log[xu]Log[xd/xu](*Phi[Sqrt[xd],Sqrt[xu],1]*))+(s+xd)Log[xd]+(s-xu)Log[xu]-4/3 ((xu-xd-1)/y)Phi[Sqrt[xd],Sqrt[xu],1]-1/3 (Log[xd]^2-Log[xu]^2)]]


(*Equation (59)*)


fCl[MC_,ml_,alpha_,mmu_,MW_,sw2_]:=(alpha^2mmu^2)/(32Pi^2MW^2sw2^2) ml^2/(MC^2-MW^2) (FCl[ml^2/MC^2]-FCl[ml^2/MW^2])
fCd[MC_,mu_,md_,alpha_,mmu_,MW_,sw2_,Qu_,Qd_,Nc_]:=(alpha^2mmu^2)/(32Pi^2MW^2sw2^2) (Nc md^2)/(MC^2-MW^2) (FCd[mu^2/MC^2,md^2/MC^2,Qu,Qd]-FCd[mu^2/MW^2,md^2/MW^2,Qu,Qd])
fCu[MC_,mu_,md_,alpha_,mmu_,MW_,sw2_,Qu_,Qd_,Nc_]:=(alpha^2mmu^2)/(32Pi^2MW^2sw2^2) (Nc mu^2)/(MC^2-MW^2) (FCu[mu^2/MC^2,md^2/MC^2,Qu+2,Qd+2]-FCu[mu^2/MW^2,md^2/MW^2,Qu+2,Qd+2])


(*Equation (58)*)


FCharged[MC_,mf_,alpha_,mmu_,MW_,sw2_,Qf_,YfA_,YmuA_]:=fCl[MC,mf[[1]],alpha,mmu,MW,sw2]YfA[[1]] YmuA+fCd[MC,mf[[3]],mf[[2]],alpha,mmu,MW,sw2,Qf[[3]],Qf[[2]],3]YfA[[2]] YmuA+fCu[MC,mf[[3]],mf[[2]],alpha,mmu,MW,sw2,Qf[[3]],Qf[[2]],3]YfA[[3]] YmuA


Fermionic2Loop[MPhi_,MhSM_,mf_,alpha_,mmu_,MW_,MZ_,sw2_,Qf_,gfv_,YfS_,YmuS_,YfhSM_]:=FNeutral[MPhi[[1]],MPhi[[2]],MPhi[[3]],mf,alpha,mmu,MW,MZ,sw2,Qf,gfv,YfS,YmuS]+
FCharged[MPhi[[4]],mf,alpha,mmu,MW,sw2,Qf,I*YfS[[3,;;]],I*YmuS[[3]]]-
Sum[( 
fgammaphi[MhSM,mf[[j]],alpha,mmu,MW,sw2,Qf[[j]],{1,3,3}[[j]]]+fZphi[MhSM,mf[[j]],alpha,mmu,MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]]
)(*YfhSM\[LeftDoubleBracket]j\[RightDoubleBracket] YfhSM\[LeftDoubleBracket]1\[RightDoubleBracket]*),{j,1,3}]
