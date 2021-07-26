#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_2loop.h"
#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/GeneralTHDM.h"
#include "gm2calc/SM.h"

#include <stdio.h>

int main()
{
   GeneralTHDM_physical_basis basis;
   basis.yukawa_scheme = GeneralTHDM_type_2;
   basis.mh = 125;
   basis.mH = 400;
   basis.mA = 420;
   basis.mHp = 440;
   basis.sin_beta_minus_alpha = 0.999;
   basis.lambda6 = 0;
   basis.lambda7 = 0;
   basis.tan_beta = 3;
   basis.m122 = 40000;

   gm2calc_SM sm;
   gm2calc_sm_set_to_default(&sm);
   sm.alpha_em_mz = 1.0/128.94579;
   sm.mu[2] = 173.34;
   sm.mu[1] = 1.28;
   sm.md[2] = 4.18;
   sm.ml[2] = 1.77684;

   GeneralTHDM* model = 0;
   gm2calc_error error = gm2calc_generalthdm_new_with_physical_basis(&model, &basis, &sm);

   if (error == gm2calc_NoError) {
      const double amu = gm2calc_generalthdm_calculate_amu_1loop(model)
                       + gm2calc_generalthdm_calculate_amu_2loop(model);

      const double delta_amu =
         gm2calc_generalthdm_calculate_uncertainty_amu_2loop(model);

      printf("amu = %g +- %g\n", amu, delta_amu);
   } else {
      printf("Error: %s\n", gm2calc_error_str(error));
   }

   gm2calc_generalthdm_free(model);

   return 0;
}
