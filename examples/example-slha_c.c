#include "gm2_1loop.h"
#include "gm2_2loop.h"
#include "MSSMNoFV_onshell.h"
#include <stdio.h>

void setup(MSSMNoFV_onshell* model) {
   /* fill pole masses */
   gm2calc_mssmnofv_set_MSvmL_pole(model, 5.18860573e+02);    /* 1L */
   gm2calc_mssmnofv_set_MSm_pole(model, 0, 5.05095249e+02);   /* 1L */
   gm2calc_mssmnofv_set_MSm_pole(model, 1, 5.25187016e+02);   /* 1L */
   gm2calc_mssmnofv_set_MChi_pole(model, 0, 2.01611468e+02);  /* 1L */
   gm2calc_mssmnofv_set_MChi_pole(model, 1, 4.10040273e+02);  /* 1L */
   gm2calc_mssmnofv_set_MChi_pole(model, 2, -5.16529941e+02); /* 1L */
   gm2calc_mssmnofv_set_MChi_pole(model, 3, 5.45628749e+02);  /* 1L */
   gm2calc_mssmnofv_set_MCha_pole(model, 0, 4.09989890e+02);  /* 1L */
   gm2calc_mssmnofv_set_MCha_pole(model, 1, 5.46057190e+02);  /* 1L */
   gm2calc_mssmnofv_set_MAh_pole(model, 1.50000000e+03);      /* 2L */

   /* fill DR-bar parameters */
   gm2calc_mssmnofv_set_TB(model, 40);                        /* 1L */
   gm2calc_mssmnofv_set_Mu(model, 500);                       /* initial guess */
   gm2calc_mssmnofv_set_MassB(model, 200);                    /* initial guess */
   gm2calc_mssmnofv_set_MassWB(model, 400);                   /* initial guess */
   gm2calc_mssmnofv_set_MassG(model, 2000);                   /* 2L */
   gm2calc_mssmnofv_set_ml2(model, 0, 0, 500 * 500);          /* 2L */
   gm2calc_mssmnofv_set_ml2(model, 1, 1, 500 * 500);          /* irrelevant */
   gm2calc_mssmnofv_set_ml2(model, 2, 2, 500 * 500);          /* 2L */
   gm2calc_mssmnofv_set_me2(model, 0, 0, 500 * 500);          /* 2L */
   gm2calc_mssmnofv_set_me2(model, 1, 1, 500 * 500);          /* initial guess */
   gm2calc_mssmnofv_set_me2(model, 2, 2, 500 * 500);          /* 2L */
   for (unsigned i = 0; i < 3; i++) {
      gm2calc_mssmnofv_set_mq2(model, i, i, 7000 * 7000);     /* 2L */
      gm2calc_mssmnofv_set_md2(model, i, i, 7000 * 7000);     /* 2L */
      gm2calc_mssmnofv_set_mu2(model, i, i, 7000 * 7000);     /* 2L */
   }
   gm2calc_mssmnofv_set_Au(model, 2, 2, 0);                   /* 2L */
   gm2calc_mssmnofv_set_Ad(model, 2, 2, 0);                   /* 2L */
   gm2calc_mssmnofv_set_Ae(model, 1, 1, 0);                   /* 1L */
   gm2calc_mssmnofv_set_Ae(model, 2, 2, 0);                   /* 2L */
   gm2calc_mssmnofv_set_scale(model, 1000);                   /* 2L */

   /* convert DR-bar parameters to on-shell */
   gm2calc_mssmnofv_convert_to_onshell(model);
}

int main() {
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();

   setup(model);

   const double amu =
      + gm2calc_calculate_amu_1loop(model)
      + gm2calc_calculate_amu_2loop(model);

   printf("amu = %e\n", amu);

   gm2calc_mssmnofv_free(model);

   return 0;
}
