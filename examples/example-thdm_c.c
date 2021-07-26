#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_2loop.h"
#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/GeneralTHDM.h"
#include "gm2calc/SM.h"

#include <stdio.h>

int main()
{
   GeneralTHDM_physical_basis* basis = gm2calc_generalthdm_physical_basis_new();
   gm2calc_generalthdm_physical_basis_set_yukawa_scheme(basis, GeneralTHDM_type_2);
   gm2calc_generalthdm_physical_basis_set_mh(basis, 125);
   gm2calc_generalthdm_physical_basis_set_mH(basis, 400);
   gm2calc_generalthdm_physical_basis_set_mA(basis, 420);
   gm2calc_generalthdm_physical_basis_set_mHp(basis, 440);
   gm2calc_generalthdm_physical_basis_set_sin_beta_minus_alpha(basis, 0.999);
   gm2calc_generalthdm_physical_basis_set_lambda_6(basis, 0);
   gm2calc_generalthdm_physical_basis_set_lambda_7(basis, 0);
   gm2calc_generalthdm_physical_basis_set_tan_beta(basis, 3);
   gm2calc_generalthdm_physical_basis_set_m122(basis, 40000);

   gm2calc_SM* sm = gm2calc_sm_new();
   gm2calc_sm_set_alpha_em_mz(sm, 1.0/128.94579);
   gm2calc_sm_set_mu(sm, 2, 173.34);
   gm2calc_sm_set_mu(sm, 1, 1.28);
   gm2calc_sm_set_md(sm, 2, 4.18);
   gm2calc_sm_set_ml(sm, 2, 1.77684);

   GeneralTHDM* model = 0;
   gm2calc_error error = gm2calc_generalthdm_new_with_physical_basis(&model, basis, sm);

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
   gm2calc_sm_free(sm);
   gm2calc_generalthdm_physical_basis_free(basis);

   return 0;
}
