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
mf = {mmu, 2, 3}; (* ml, md, mu *)
Qf = { -1, -1/3, 2/3 };
YmuA = 3;
YfA = { YmuA, 4, 5 };

amuC = 3 FCharged[MC, mf, alpha, mmu, MW, sw2, Qf, YfA, YmuA]

Print["amuC 10^10 = ", N[amuC 10^10, 17]]

*/

TEST_CASE("2-loop_fermionic_charged")
{
   gm2calc::general_thdm::THDM_F_parameters pars;
   pars.alpha = 1./137;
   pars.mw = 80.0;
   pars.mz = 90.0;
   pars.mHp = 200.0;
   pars.mu << 0.3, 3.0, 30.0;
   pars.md << 0.2, 2.0, 20.0;
   pars.ml << 0.1, 1.0, 10.0;
   pars.yuS << 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0;
   pars.ydS << 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0;
   pars.ylS << 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0;

   const auto amu = gm2calc::general_thdm::amu2L_F_charged(pars);

   CHECK_CLOSE(1e10*amu, 33.045119891312677, 1e-12);
}

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
mf = {mmu, 2, 3}; (* ml, md, mu *)
Qf = { -1, -1/3, 2/3 }; (* l, d, u *)
YmuA = 3;
YfA = { YmuA, 4, 5 };
MhSM = 125;
Mphi = { 100, 400, 300 }; (* h, H, A *)
T3 = { -1, -1, 1 }; (* l, d, u *)
gfv = T3/2 - Qf sw2;
Yfphi = { YfA, YfA, YfA };
Ymuphi = { YmuA, YmuA, YmuA };

amuN = 3 (
    FNeutral[Mphi,mf,alpha,mmu,MW,MZ,sw2,Qf,gfv,Yfphi,Ymuphi] -
    Sum[fgammaphi[MhSM,mf[[j]],alpha,mmu,MW,sw2,Qf[[j]],{1,3,3}[[j]]]
        + fZphi[MhSM,mf[[j]],alpha,mmu,MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]],
        {j,1,3}
    ])

Print["amuN 10^10 = ", N[10^10 amuN, 17]]

*/

TEST_CASE("2-loop_fermionic_neutral")
{
   gm2calc::general_thdm::THDM_F_parameters pars;
   pars.alpha = 1./137;
   pars.mw = 80.0;
   pars.mz = 90.0;
   pars.mh << 100.0, 400.0;
   pars.mhSM = 125.0;
   pars.mA = 300.0;
   pars.mHp = 200.0;
   pars.mu << 0.3, 3.0, 30.0;
   pars.md << 0.2, 2.0, 20.0;
   pars.ml << 0.1, 1.0, 10.0;
   pars.yuS << 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0;
   pars.ydS << 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0;
   pars.ylS << 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0;

   const auto amu = gm2calc::general_thdm::amu2L_F_neutral(pars);

   CHECK_CLOSE(1e10*amu, -100.23855413383224, 1e-12);
}
