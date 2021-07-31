#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_2loop.h"
#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/THDM.h"
#include "gm2calc/SM.h"

#include <stdio.h>

int main()
{
   THDM_mass_basis basis;
   basis.yukawa_type = THDM_type_2;
   basis.mh = 125;
   basis.mH = 400;
   basis.mA = 420;
   basis.mHp = 440;
   basis.sin_beta_minus_alpha = 0.999;
   basis.lambda_6 = 0;
   basis.lambda_7 = 0;
   basis.tan_beta = 3;
   basis.m122 = 40000;
   basis.zeta_u = 0;
   basis.zeta_d = 0;
   basis.zeta_l = 0;
   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         basis.Xu_real[i][k] = 0;
         basis.Xu_imag[i][k] = 0;
         basis.Xd_real[i][k] = 0;
         basis.Xd_imag[i][k] = 0;
         basis.Xl_real[i][k] = 0;
         basis.Xl_imag[i][k] = 0;
      }
   }

   SM sm;
   gm2calc_sm_set_to_default(&sm);
   sm.alpha_em_mz = 1.0/128.94579;
   sm.mu[2] = 173.34;
   sm.mu[1] = 1.28;
   sm.md[2] = 4.18;
   sm.ml[2] = 1.77684;

   THDM_config config;
   config.running_couplings = 0;

   THDM* model = 0;
   gm2calc_error error = gm2calc_thdm_new_with_mass_basis(&model, &basis, &sm, &config);

   if (error == gm2calc_NoError) {
      const double amu = gm2calc_thdm_calculate_amu_1loop(model)
                       + gm2calc_thdm_calculate_amu_2loop(model);

      const double delta_amu =
         gm2calc_thdm_calculate_uncertainty_amu_2loop(model);

      printf("amu = %g +- %g\n", amu, delta_amu);
   } else {
      printf("Error: %s\n", gm2calc_error_str(error));
   }

   gm2calc_thdm_free(model);

   return 0;
}
