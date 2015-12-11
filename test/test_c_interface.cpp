#include "test.hpp"

#include "gm2_1loop.h"
#include "gm2_2loop.h"
#include "MSSMNoFV_onshell.h"

#include "gm2_1loop.hpp"
#include "gm2_2loop.hpp"
#include "MSSMNoFV_onshell.hpp"

using namespace gm2calc_test;

void setup(MSSMNoFV_onshell* model)
{
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
   gm2calc_mssmnofv_set_MAh_pole(model,1500); /* 2L */
   gm2calc_mssmnofv_set_scale(model,454.7); /* 2L */

   for (unsigned i = 0; i < 3; i++) {
      gm2calc_mssmnofv_set_mq2(model,i,i,500*500); /* 2L */
      gm2calc_mssmnofv_set_ml2(model,i,i,500*500); /* 1L(smuon)/2L */
      gm2calc_mssmnofv_set_md2(model,i,i,500*500); /* 2L */
      gm2calc_mssmnofv_set_mu2(model,i,i,500*500); /* 2L */
      gm2calc_mssmnofv_set_me2(model,i,i,500*500); /* 1L(smuon)/2L */
   }

   /* calculate mass spectrum */
   gm2calc_mssmnofv_calculate_masses(model);
}

void test_parameters(MSSMNoFV_onshell* model)
{
   const gm2calc::MSSMNoFV_onshell mcpp(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));

#define COMPARE_0(a)                                                    \
   CHECK_EQUAL(gm2calc_mssmnofv_get_ ## a(model), mcpp.get_ ## a())
#define COMPARE_1(a,i)                                                  \
   CHECK_EQUAL(gm2calc_mssmnofv_get_ ## a(model,i), mcpp.get_ ## a(i))
#define COMPARE_2(a,i,k)                                                \
   CHECK_EQUAL(gm2calc_mssmnofv_get_ ## a(model,i,k), mcpp.get_ ## a(i,k))

   for (unsigned i = 0; i < 3; i++) {
      for (unsigned k = 0; k < 3; k++) {
         COMPARE_2(Ae,i,k);
         COMPARE_2(Ad,i,k);
         COMPARE_2(Au,i,k);
         COMPARE_2(mq2,i,k);
         COMPARE_2(md2,i,k);
         COMPARE_2(mu2,i,k);
         COMPARE_2(ml2,i,k);
         COMPARE_2(me2,i,k);
         COMPARE_2(Ye,i,k);
         COMPARE_2(Yd,i,k);
         COMPARE_2(Yu,i,k);
      }
   }

   COMPARE_0(EL);
   COMPARE_0(EL0);
   COMPARE_0(gY);
   COMPARE_0(g1);
   COMPARE_0(g2);
   COMPARE_0(g3);
   COMPARE_0(TB);
   COMPARE_0(MassB);
   COMPARE_0(MassWB);
   COMPARE_0(MassG);
   COMPARE_0(Mu);
   COMPARE_0(vev);
   COMPARE_0(MW);
   COMPARE_0(MZ);
   COMPARE_0(ME);
   COMPARE_0(MM);
   COMPARE_0(ML);
   COMPARE_0(MU);
   COMPARE_0(MC);
   COMPARE_0(MT);
   COMPARE_0(MD);
   COMPARE_0(MS);
   COMPARE_0(MB);
   COMPARE_0(MBMB);
   COMPARE_1(MCha,0);
   COMPARE_1(MCha,1);
   COMPARE_1(MChi,0);
   COMPARE_1(MChi,1);
   COMPARE_1(MChi,2);
   COMPARE_1(MChi,3);
   COMPARE_1(MSm,0);
   COMPARE_1(MSm,1);
   COMPARE_0(MSvmL);

   CHECK_EQUAL(gm2calc_mssmnofv_get_MAh(model), mcpp.get_MAh(1));

#undef COMPARE_0
#undef COMPARE_1
#undef COMPARE_2
}

void test_1_loop(MSSMNoFV_onshell* model)
{
   const gm2calc::MSSMNoFV_onshell mcpp(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));

   CHECK_EQUAL(gm2calc_amuChi0(model), gm2calc::amuChi0(mcpp));
   CHECK_EQUAL(gm2calc_amuChipm(model), gm2calc::amuChipm(mcpp));

   CHECK_EQUAL(gm2calc_calculate_amu_1loop(model), gm2calc::calculate_amu_1loop(mcpp));
   CHECK_EQUAL(gm2calc_calculate_amu_1loop_non_tan_beta_resummed(model),
               gm2calc::calculate_amu_1loop_non_tan_beta_resummed(mcpp));
}

void test_2_loop(MSSMNoFV_onshell* model)
{
   const gm2calc::MSSMNoFV_onshell mcpp(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));

   // fermion/sfermion 2L corrections
   CHECK_EQUAL(gm2calc_amu2LFSfapprox(model), gm2calc::amu2LFSfapprox(mcpp));
   CHECK_EQUAL(gm2calc_amu2LFSfapprox_non_tan_beta_resummed(model),
               gm2calc::amu2LFSfapprox_non_tan_beta_resummed(mcpp));

   // photonic 2L corrections
   CHECK_EQUAL(gm2calc_amuChi0Photonic(model), gm2calc::amuChi0Photonic(mcpp));
   CHECK_EQUAL(gm2calc_amuChipmPhotonic(model), gm2calc::amuChipmPhotonic(mcpp));

   // 2L(a) diagrams
   CHECK_EQUAL(gm2calc_amu2LaSferm(model), gm2calc::amu2LaSferm(mcpp));
   CHECK_EQUAL(gm2calc_amu2LaCha(model), gm2calc::amu2LaCha(mcpp));

   CHECK_EQUAL(gm2calc_calculate_amu_2loop(model), gm2calc::calculate_amu_2loop(mcpp));
   CHECK_EQUAL(gm2calc_calculate_amu_2loop_non_tan_beta_resummed(model),
               gm2calc::calculate_amu_2loop_non_tan_beta_resummed(mcpp));
}

int main()
{
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();

   setup(model);

   print_mssmnofv(model);

   printf("\n");
   printf("==============================\n");
   printf("testing parameter getters\n");
   printf("==============================\n");
   printf("\n");

   test_parameters(model);

   printf("\n");
   printf("==============================\n");
   printf("testing amu 1-loop calculation\n");
   printf("==============================\n");
   printf("\n");

   test_1_loop(model);

   printf("\n");
   printf("==============================\n");
   printf("testing amu 2-loop calculation\n");
   printf("==============================\n");
   printf("\n");

   test_2_loop(model);

   printf("\n");
   printf("==============================\n");
   printf("Test results: %s\n", (g_error ? "FAIL" : "OK"));
   printf("==============================\n");

   gm2calc_mssmnofv_free(model);

   return g_error;
}
