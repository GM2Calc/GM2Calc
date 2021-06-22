#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "GeneralTHDM/gm2_2loop_helpers.hpp"

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


TEST_CASE("2-loop_fermionic_charged")
{
   gm2calc::general_thdm::THDM_F_parameters pars;
   const auto amu = gm2calc::general_thdm::amu2L_F_charged(pars);

   // CHECK_CLOSE(amu, 0.0, 1e-15);
}
