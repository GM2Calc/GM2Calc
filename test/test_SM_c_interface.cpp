#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/SM.hpp"
#include "gm2calc/SM.h"


#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


TEST_CASE("default")
{
   const gm2calc::SM sm;

   gm2calc_SM sm_c;
   gm2calc_sm_set_to_default(&sm_c);

   CHECK_EQ(sm_c.alpha_em_0, sm.get_alpha_em_0());
   CHECK_EQ(sm_c.alpha_em_mz, sm.get_alpha_em_mz());
   CHECK_EQ(sm_c.alpha_s_mz, sm.get_alpha_s_mz());
   CHECK_EQ(sm_c.mh, sm.get_mh());
   CHECK_EQ(sm_c.mw, sm.get_mw());
   CHECK_EQ(sm_c.mz, sm.get_mz());

   for (int i = 0; i < 3; i++) {
      CHECK_EQ(sm_c.mu[i], sm.get_mu(i));
      CHECK_EQ(sm_c.md[i], sm.get_md(i));
      CHECK_EQ(sm_c.mv[i], sm.get_mv(i));
      CHECK_EQ(sm_c.ml[i], sm.get_ml(i));
   }

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         CHECK_EQ(sm_c.ckm_real[i][k], std::real(sm.get_ckm()(i, k)));
         CHECK_EQ(sm_c.ckm_imag[i][k], std::imag(sm.get_ckm()(i, k)));
      }
   }
}
