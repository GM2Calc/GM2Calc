#include "gm2_1loop.h"
#include "gm2_2loop.h"
#include "MSSMNoFV_onshell.h"

#include "gm2_1loop.hpp"
#include "gm2_2loop.hpp"
#include "MSSMNoFV_onshell.hpp"

#include <iostream>
#include <string>

static int g_error = 0;
static bool verbose = false;

#define VERBOSE_MSG(...)                        \
   do {                                         \
   if (verbose) {                               \
      std::cout << __VA_ARGS__ << std::endl;    \
   }                                            \
   } while (0)

#define CHECK_EQUAL(a, b) \
   check_equal((a), (b), #a, #b)

void check_equal(double a, double b, const std::string& astr, const std::string& bstr)
{
   if (a == b) {
      VERBOSE_MSG("test passed: " << astr << " == " << bstr);
   } else {
      std::cout << "test failed: " << astr << " == " << bstr
                << ": " << a << " != " << b << std::endl;
      g_error = 1;
   }
}

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
   gm2calc_mssmnofv_set_MA0(model,1500);    /* 2L */
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

int main()
{
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();

   setup(model);

   print_mssmnofv(model);

   printf("\n");
   printf("==============================\n");
   printf("testing amu 1-loop calculation\n");
   printf("==============================\n");
   printf("\n");

   test_1_loop(model);

   printf("\n");
   printf("==============================\n");
   printf("Test results: %s\n", (g_error ? "FAIL" : "OK"));
   printf("==============================\n");

   gm2calc_mssmnofv_free(model);

   return g_error;
}
