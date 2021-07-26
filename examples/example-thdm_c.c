#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_2loop.h"
#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/GeneralTHDM.h"
#include "gm2calc/SM.h"

#include <stdio.h>

int main()
{
   GeneralTHDM_physical_basis* basis = gm2calc_generalthdm_physical_basis_new();
   gm2calc_SM* sm = gm2calc_generalthdm_sm_new();

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
   gm2calc_generalthdm_sm_free(sm);
   gm2calc_generalthdm_physical_basis_free(basis);

   return 0;
}
