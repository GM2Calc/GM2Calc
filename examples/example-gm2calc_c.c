#include "gm2_1loop.h"
#include "gm2_2loop.h"
#include "MSSMNoFV_onshell.h"
#include <stdio.h>
#include <stdlib.h>

void setup(MSSMNoFV_onshell* model) {
   /* fill DR-bar parameters */
   gm2calc_mssmnofv_set_TB(model, 10);      /* 1L */
   gm2calc_mssmnofv_set_Ae(model,1,1,0);    /* 1L */

   /* fill on-shell parameters */
   gm2calc_mssmnofv_set_Mu(model,350);      /* 1L */
   gm2calc_mssmnofv_set_MassB(model,150);   /* 1L */
   gm2calc_mssmnofv_set_MassWB(model,300);  /* 1L */
   gm2calc_mssmnofv_set_MassG(model,1000);  /* 2L */
   gm2calc_mssmnofv_set_Au(model,2,2,0);    /* 2L */
   gm2calc_mssmnofv_set_Ad(model,2,2,0);    /* 2L */
   gm2calc_mssmnofv_set_Ae(model,2,2,0);    /* 2L */
   gm2calc_mssmnofv_set_MAh_pole(model,1500);    /* 2L */
   gm2calc_mssmnofv_set_scale(model,454.7); /* 2L */

   for (unsigned i = 0; i < 3; i++) {
      gm2calc_mssmnofv_set_mq2(model,i,i,500*500); /* 2L */
      gm2calc_mssmnofv_set_ml2(model,i,i,500*500); /* 1L(smuon)/2L */
      gm2calc_mssmnofv_set_md2(model,i,i,500*500); /* 2L */
      gm2calc_mssmnofv_set_mu2(model,i,i,500*500); /* 2L */
      gm2calc_mssmnofv_set_me2(model,i,i,500*500); /* 1L(smuon)/2L */
   }

   /* calculate mass spectrum */
   const enum EError error = gm2calc_mssmnofv_calculate_masses(model);

   if (error != NoError) {
      printf("Error: %s\n", gm2calc_error_str(error));
      abort();
   }
}

int main() {
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();

   setup(model);

   const double amu =
      + gm2calc_calculate_amu_1loop(model)
      + gm2calc_calculate_amu_2loop(model);

   printf("amu = %e\n", amu);

   /* destroy model to prevent resource leak */
   gm2calc_mssmnofv_free(model);

   return 0;
}
