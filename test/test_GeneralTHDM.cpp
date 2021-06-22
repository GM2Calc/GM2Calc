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

Print[N[amuC, 17]]

*/

TEST_CASE("2-loop_fermionic_charged")
{
   gm2calc::general_thdm::THDM_F_parameters pars;
   pars.alpha = 1./137;
   pars.mw = 80.0;
   pars.mz = 90.0;
   pars.mHp = 200.0;
   pars.mu << 3.0, 3.0, 3.0;
   pars.md << 2.0, 2.0, 2.0;
   pars.ml << 1.0, 1.0, 1.0;
   pars.yuS << 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0;
   pars.ydS << 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0;
   pars.ylS << 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0;

   const auto amu = gm2calc::general_thdm::amu2L_F_charged(pars);

   CHECK_CLOSE(1e10*amu, 3.0942542015032706, 1e-12);
}
