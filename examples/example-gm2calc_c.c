/* #include "gm2_1loop.h" */
/* #include "gm2_2loop.h" */
#include "MSSMNoFV_onshell.h"
#include <stdio.h>

void setup(C_MSSMNoFV_onshell* model) {
   /* fill DR-bar parameters */
   c_mssm_onshell_set_TB(model, 10);      /* 1L */
   c_mssm_onshell_set_Ae(model,1,1,0);    /* 1L */

   /* fill on-shell parameters */
   c_mssm_onshell_set_Mu(model,350);      /* 1L */
   c_mssm_onshell_set_MassB(model,150);   /* 1L */
   c_mssm_onshell_set_MassWB(model,300);  /* 1L */
   c_mssm_onshell_set_MassG(model,1000);  /* 2L */
   c_mssm_onshell_set_mq2(model,500*500); /* 2L */
   c_mssm_onshell_set_ml2(model,500*500); /* 1L(smuon)/2L */
   c_mssm_onshell_set_md2(model,500*500); /* 2L */
   c_mssm_onshell_set_mu2(model,500*500); /* 2L */
   c_mssm_onshell_set_me2(model,500*500); /* 1L(smuon)/2L */
   c_mssm_onshell_set_Au(model,2,2,0);    /* 2L */
   c_mssm_onshell_set_Ad(model,2,2,0);    /* 2L */
   c_mssm_onshell_set_Ae(model,2,2,0);    /* 2L */
   c_mssm_onshell_set_MA0(model,1500);    /* 2L */
   c_mssm_onshell_set_scale(model,454.7); /* 2L */

   /* calculate mass spectrum */
   c_mssm_onshell_calculate_masses(model);
}

int main() {
   double amu;
   C_MSSMNoFV_onshell* model = c_mssm_onshell_new();

   setup(model);

   amu = c_calculate_amu_1loop(model) + c_calculate_amu_2loop(model);

   printf("amu = %e\n", amu);

   c_mssm_onshell_free(model);

   return 0;
}
