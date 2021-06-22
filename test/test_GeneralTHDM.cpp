#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "GeneralTHDM/gm2_2loop_helpers.hpp"

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


/* Generated with:

<< ../math/GeneralTHDMBarrZee.m;

MW = 80;
MZ = 90;
cw = MW/MZ;
cw2 = cw^2;
sw2 = 1 - cw2;
MC = 200;
alpha = 1/137;
mmu = 1;
mf1 = {mmu, 2, 3}/10; (* me, md, mu *)
mf2 = {mmu, 2, 3};    (* mm, ms, mc *)
mf3 = {mmu, 2, 3}*10; (* mtau, mb, mt *)
Qf = { -1, -1/3, 2/3 }; (* l, d, u *)
YmuA = 3;
Yf1A = { YmuA, 4, 5 }/10;
Yf2A = { YmuA, 4, 5 };
Yf3A = { YmuA, 4, 5 }*10;

amuC = (
    FCharged[MC, mf1, alpha, mmu, MW, sw2, Qf, Yf1A, YmuA] +
    FCharged[MC, mf2, alpha, mmu, MW, sw2, Qf, Yf2A, YmuA] +
    FCharged[MC, mf3, alpha, mmu, MW, sw2, Qf, Yf3A, YmuA]
)

Print["amuC 10^10 = ", N[10^10 amuC, 17]]

*/

TEST_CASE("2-loop_fermionic_charged")
{
   gm2calc::general_thdm::THDM_F_parameters pars;
   pars.alpha = 1./137;
   pars.mm = 1.0;
   pars.mw = 80.0;
   pars.mz = 90.0;
   pars.mHp = 200.0;
   pars.mu << 0.3, 3.0, 30.0;
   pars.md << 0.2, 2.0, 20.0;
   pars.ml << 0.1, 1.0, 10.0;
   pars.yuS << 0.5, 0.5, 0.5, 5.0, 5.0, 5.0, 50.0, 50.0, 50.0;
   pars.ydS << 0.4, 0.4, 0.4, 4.0, 4.0, 4.0, 40.0, 40.0, 40.0;
   pars.ylS << 0.3, 0.3, 0.3, 3.0, 3.0, 3.0, 30.0, 30.0, 30.0;

   const auto amu = gm2calc::general_thdm::amu2L_F_charged(pars);

   CHECK_CLOSE(1e10*amu, 320.98648556040420, 1e-12);
}

/* Generated with:

MW = 80;
MZ = 90;
cw = MW/MZ;
cw2 = cw^2;
sw2 = 1 - cw2;
MC = 200;
alpha = 1/137;
mmu = 1;
mf1 = {mmu, 2, 3}/10; (* me, md, mu *)
mf2 = {mmu, 2, 3};    (* mm, ms, mc *)
mf3 = {mmu, 2, 3}*10; (* mtau, mb, mt *)
Qf = { -1, -1/3, 2/3 }; (* l, d, u *)
YmuA = 3;
Yf1A = { YmuA, 4, 5 }/10;
Yf2A = { YmuA, 4, 5 };
Yf3A = { YmuA, 4, 5 }*10;
MhSM = 125;
Mphi = { 100, 400, 300 }; (* h, H, A *)
T3 = { -1, -1, 1 }; (* l, d, u *)
gfv = T3/2 - Qf sw2;
Yf1phi = { Yf1A, Yf1A, Yf1A };
Yf2phi = { Yf2A, Yf2A, Yf2A };
Yf3phi = { Yf3A, Yf3A, Yf3A };
Ymuphi = { YmuA, YmuA, YmuA };

amuN = (
    FNeutral[Mphi,mf1,alpha,mmu,MW,MZ,sw2,Qf,gfv,Yf1phi,Ymuphi] -
    Sum[fgammaphi[MhSM,mf1[[j]],alpha,mmu,MW,sw2,Qf[[j]],{1,3,3}[[j]]]
        + fZphi[MhSM,mf1[[j]],alpha,mmu,MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]],
        {j,1,3}] +
    FNeutral[Mphi,mf2,alpha,mmu,MW,MZ,sw2,Qf,gfv,Yf2phi,Ymuphi] -
    Sum[fgammaphi[MhSM,mf2[[j]],alpha,mmu,MW,sw2,Qf[[j]],{1,3,3}[[j]]]
        + fZphi[MhSM,mf2[[j]],alpha,mmu,MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]],
        {j,1,3}] +
    FNeutral[Mphi,mf3,alpha,mmu,MW,MZ,sw2,Qf,gfv,Yf3phi,Ymuphi] -
    Sum[fgammaphi[MhSM,mf3[[j]],alpha,mmu,MW,sw2,Qf[[j]],{1,3,3}[[j]]]
        + fZphi[MhSM,mf3[[j]],alpha,mmu,MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]],
        {j,1,3}]
)

Print["amuN 10^10 = ", N[10^10 amuN, 17]]

*/

TEST_CASE("2-loop_fermionic_neutral")
{
   gm2calc::general_thdm::THDM_F_parameters pars;
   pars.alpha = 1./137;
   pars.mm = 1.0;
   pars.mw = 80.0;
   pars.mz = 90.0;
   pars.mh << 100.0, 400.0;
   pars.mhSM = 125.0;
   pars.mA = 300.0;
   pars.mHp = 200.0;
   pars.mu << 0.3, 3.0, 30.0;
   pars.md << 0.2, 2.0, 20.0;
   pars.ml << 0.1, 1.0, 10.0;
   pars.yuS << 0.5, 0.5, 0.5, 5.0, 5.0, 5.0, 50.0, 50.0, 50.0;
   pars.ydS << 0.4, 0.4, 0.4, 4.0, 4.0, 4.0, 40.0, 40.0, 40.0;
   pars.ylS << 0.3, 0.3, 0.3, 3.0, 3.0, 3.0, 30.0, 30.0, 30.0;

   const auto amu = gm2calc::general_thdm::amu2L_F_neutral(pars);

   CHECK_CLOSE(1e10*amu, -1023.0630847112763, 1e-12);
}


TEST_CASE("2-loop_bosonic")
{
   gm2calc::general_thdm::THDM_B_parameters pars;
   pars.alpha = 1./137;
   pars.mm = 1.0;
   pars.mw = 80.0;
   pars.mz = 90.0;
   pars.mh << 100.0, 400.0;
   pars.mhSM = 125.0;
   pars.mA = 300.0;
   pars.mHp = 210.0;
   pars.tb = 30.0;
   pars.zetal = 2.0;
   pars.eta = 3.0;
   pars.lambda5 = 4.0;

   const auto amu = gm2calc::general_thdm::amu2L_B(pars);

   // CHECK_CLOSE(1e10*amu, 0.0, 1e-12);
}
