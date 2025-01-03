#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"

#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_2loop.h"
#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/MSSMNoFV_onshell.h"

#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/MSSMNoFV_onshell.hpp"

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


#define CHECK_EQUAL(a,b)                                \
   CHECK_CLOSE(a,b,std::numeric_limits<double>::epsilon())


class Cleanup_on_destruction {
public:
   Cleanup_on_destruction(MSSMNoFV_onshell* model_) : model(model_) {}
   ~Cleanup_on_destruction() { gm2calc_mssmnofv_free(model); }
private:
   MSSMNoFV_onshell* model = nullptr;
};


void setup_gm2calc_scheme(MSSMNoFV_onshell* model)
{
   /* fill SM parameters */
   gm2calc_mssmnofv_set_alpha_MZ(model, 1./127);
   gm2calc_mssmnofv_set_alpha_thompson(model, 1./137);
   gm2calc_mssmnofv_set_g3(model, 0.118);
   gm2calc_mssmnofv_set_MZ_pole(model, 91);
   gm2calc_mssmnofv_set_MW_pole(model, 83);
   gm2calc_mssmnofv_set_MT_pole(model, 172);
   gm2calc_mssmnofv_set_MB_running(model, 3);
   gm2calc_mssmnofv_set_ML_pole(model, 1.7);
   gm2calc_mssmnofv_set_MM_pole(model, 0.1);

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


gm2calc::MSSMNoFV_onshell setup_gm2calc_scheme()
{
   gm2calc::MSSMNoFV_onshell model;

   const Eigen::Matrix<double,3,3> UnitMatrix
      = Eigen::Matrix<double,3,3>::Identity();

   // fill SM parameters
   model.set_alpha_MZ(1./127);
   model.set_alpha_thompson(1./137);
   model.set_g3(0.118);
   model.get_physical().MVZ = 91;
   model.get_physical().MVWm = 83;
   model.get_physical().MFt = 172;
   model.get_physical().MFb = 3;
   model.get_physical().MFtau = 1.7;
   model.get_physical().MFm = 0.1;

   // fill DR-bar parameters
   model.set_TB(10);                      // 1L
   model.set_Ae(1,1,0);                   // 1L

   // fill on-shell parameters
   model.set_Mu(350);                     // 1L
   model.set_MassB(150);                  // 1L
   model.set_MassWB(300);                 // 1L
   model.set_MassG(1000);                 // 2L
   model.set_mq2(500 * 500 * UnitMatrix); // 2L
   model.set_ml2(500 * 500 * UnitMatrix); // 1L(smuon)/2L
   model.set_md2(500 * 500 * UnitMatrix); // 2L
   model.set_mu2(500 * 500 * UnitMatrix); // 2L
   model.set_me2(500 * 500 * UnitMatrix); // 1L(smuon)/2L
   model.set_Au(2,2,0);                   // 2L
   model.set_Ad(2,2,0);                   // 2L
   model.set_Ae(2,2,0);                   // 2L
   model.set_MA0(1500);                   // 2L
   model.set_scale(454.7);                // 2L

   // calculate mass spectrum
   model.calculate_masses();

   return model;
}


void setup_slha_scheme(MSSMNoFV_onshell* model)
{
   setup_gm2calc_scheme(model);

   gm2calc_mssmnofv_set_MSm_pole(model, 0, 300);
   gm2calc_mssmnofv_set_MSm_pole(model, 1, 400);
   gm2calc_mssmnofv_set_MSvmL_pole(model, 500);
   gm2calc_mssmnofv_set_MCha_pole(model, 0, 600);
   gm2calc_mssmnofv_set_MCha_pole(model, 1, 700);
   gm2calc_mssmnofv_set_MChi_pole(model, 0, 800);
   gm2calc_mssmnofv_set_MChi_pole(model, 1, 900);
   gm2calc_mssmnofv_set_MChi_pole(model, 2, 1000);
   gm2calc_mssmnofv_set_MChi_pole(model, 3, 1100);

   gm2calc_mssmnofv_convert_to_onshell(model);
}


gm2calc::MSSMNoFV_onshell setup_slha_scheme()
{
   auto model = setup_gm2calc_scheme();

   model.get_physical().MSm(0) = 300;
   model.get_physical().MSm(1) = 400;
   model.get_physical().MSvmL = 500;
   model.get_physical().MCha(0) = 600;
   model.get_physical().MCha(1) = 700;
   model.get_physical().MChi(0) = 800;
   model.get_physical().MChi(1) = 900;
   model.get_physical().MChi(2) = 1000;
   model.get_physical().MChi(3) = 1100;

   model.convert_to_onshell();

   return model;
}


void test_parameters(const MSSMNoFV_onshell* model, const gm2calc::MSSMNoFV_onshell& model2)
{
#define COMPARE_0(a)                                                    \
   CHECK_EQUAL(gm2calc_mssmnofv_get_ ## a(model), model2.get_ ## a())
#define COMPARE_1(a,i)                                                  \
   CHECK_EQUAL(gm2calc_mssmnofv_get_ ## a(model,i), model2.get_ ## a(i))
#define COMPARE_2(a,i,k)                                                \
   CHECK_EQUAL(gm2calc_mssmnofv_get_ ## a(model,i,k), model2.get_ ## a(i,k))

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

   for (unsigned i = 0; i < 2; i++) {
      for (unsigned k = 0; k < 2; k++) {
         CHECK_EQUAL(gm2calc_mssmnofv_get_USe(model,i,k)  , model2.get_ZE(i,k));
         CHECK_EQUAL(gm2calc_mssmnofv_get_USm(model,i,k)  , model2.get_ZM(i,k));
         CHECK_EQUAL(gm2calc_mssmnofv_get_UStau(model,i,k), model2.get_ZTau(i,k));
         CHECK_EQUAL(gm2calc_mssmnofv_get_USu(model,i,k)  , model2.get_ZU(i,k));
         CHECK_EQUAL(gm2calc_mssmnofv_get_USd(model,i,k)  , model2.get_ZD(i,k));
         CHECK_EQUAL(gm2calc_mssmnofv_get_USc(model,i,k)  , model2.get_ZC(i,k));
         CHECK_EQUAL(gm2calc_mssmnofv_get_USs(model,i,k)  , model2.get_ZS(i,k));
         CHECK_EQUAL(gm2calc_mssmnofv_get_USt(model,i,k)  , model2.get_ZT(i,k));
         CHECK_EQUAL(gm2calc_mssmnofv_get_USb(model,i,k)  , model2.get_ZB(i,k));

         double im_UM = 0.0, im_UP = 0.0;
         const double re_UM = gm2calc_mssmnofv_get_UM(model, i, k, &im_UM);
         const double re_UP = gm2calc_mssmnofv_get_UP(model, i, k, &im_UP);

         CHECK_EQUAL(re_UM, std::real(model2.get_UM(i,k)));
         CHECK_EQUAL(im_UM, std::imag(model2.get_UM(i,k)));
         CHECK_EQUAL(re_UP, std::real(model2.get_UP(i,k)));
         CHECK_EQUAL(im_UP, std::imag(model2.get_UP(i,k)));
      }
   }

   for (unsigned i = 0; i < 4; i++) {
      for (unsigned k = 0; k < 4; k++) {
         double im_ZN = 0.0;
         const double re_ZN = gm2calc_mssmnofv_get_ZN(model, i, k, &im_ZN);
         CHECK_EQUAL(re_ZN, std::real(model2.get_ZN(i,k)));
         CHECK_EQUAL(im_ZN, std::imag(model2.get_ZN(i,k)));
      }
   }

   COMPARE_0(scale);
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
   COMPARE_1(MSu,0);
   COMPARE_1(MSu,1);
   COMPARE_1(MSd,0);
   COMPARE_1(MSd,1);
   COMPARE_1(MSc,0);
   COMPARE_1(MSc,1);
   COMPARE_1(MSs,0);
   COMPARE_1(MSs,1);
   COMPARE_1(MSt,0);
   COMPARE_1(MSt,1);
   COMPARE_1(MSb,0);
   COMPARE_1(MSb,1);
   COMPARE_1(MSe,0);
   COMPARE_1(MSe,1);
   COMPARE_1(MSm,0);
   COMPARE_1(MSm,1);
   COMPARE_1(MSe,0);
   COMPARE_1(MSe,1);
   COMPARE_1(MStau,0);
   COMPARE_1(MStau,1);
   COMPARE_0(MSveL);
   COMPARE_0(MSvmL);
   COMPARE_0(MSvtL);
   COMPARE_1(Mhh,0);
   COMPARE_1(Mhh,1);

   CHECK_EQUAL(gm2calc_mssmnofv_get_MAh(model), model2.get_MAh(1));

#undef COMPARE_0
#undef COMPARE_1
#undef COMPARE_2
}


TEST_CASE("parameter_setters")
{
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   Cleanup_on_destruction cleanup(model);

   setup_gm2calc_scheme(model);
   const gm2calc::MSSMNoFV_onshell model2 = setup_gm2calc_scheme();

   test_parameters(model, model2);
}


TEST_CASE("parameter_getters")
{
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   Cleanup_on_destruction cleanup(model);

   setup_gm2calc_scheme(model);

   const gm2calc::MSSMNoFV_onshell mcpp(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));

   test_parameters(model, mcpp);
}


TEST_CASE("conversion")
{
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   Cleanup_on_destruction cleanup(model);

   setup_slha_scheme(model);
   const gm2calc::MSSMNoFV_onshell model2 = setup_slha_scheme();

   test_parameters(model, model2);
}


TEST_CASE("1_loop")
{
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   Cleanup_on_destruction cleanup(model);

   setup_gm2calc_scheme(model);

   const gm2calc::MSSMNoFV_onshell mcpp(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));

   CHECK_EQUAL(gm2calc_mssmnofv_amu1LChi0(model), gm2calc::amu1LChi0(mcpp));
   CHECK_EQUAL(gm2calc_mssmnofv_amu1LChipm(model), gm2calc::amu1LChipm(mcpp));

   CHECK_EQUAL(gm2calc_mssmnofv_calculate_amu_1loop(model), gm2calc::calculate_amu_1loop(mcpp));
   CHECK_EQUAL(gm2calc_mssmnofv_calculate_amu_1loop_non_tan_beta_resummed(model),
               gm2calc::calculate_amu_1loop_non_tan_beta_resummed(mcpp));
}


TEST_CASE("2_loop")
{
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   Cleanup_on_destruction cleanup(model);

   setup_gm2calc_scheme(model);

   const gm2calc::MSSMNoFV_onshell mcpp(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));

   // fermion/sfermion 2L corrections
   CHECK_EQUAL(gm2calc_mssmnofv_amu2LFSfapprox(model), gm2calc::amu2LFSfapprox(mcpp));
   CHECK_EQUAL(gm2calc_mssmnofv_amu2LFSfapprox_non_tan_beta_resummed(model),
               gm2calc::amu2LFSfapprox_non_tan_beta_resummed(mcpp));

   // photonic 2L corrections
   CHECK_EQUAL(gm2calc_mssmnofv_amu2LChi0Photonic(model), gm2calc::amu2LChi0Photonic(mcpp));
   CHECK_EQUAL(gm2calc_mssmnofv_amu2LChipmPhotonic(model), gm2calc::amu2LChipmPhotonic(mcpp));

   // 2L(a) diagrams
   CHECK_EQUAL(gm2calc_mssmnofv_amu2LaSferm(model), gm2calc::amu2LaSferm(mcpp));
   CHECK_EQUAL(gm2calc_mssmnofv_amu2LaCha(model), gm2calc::amu2LaCha(mcpp));

   CHECK_EQUAL(gm2calc_mssmnofv_calculate_amu_2loop(model), gm2calc::calculate_amu_2loop(mcpp));
   CHECK_EQUAL(gm2calc_mssmnofv_calculate_amu_2loop_non_tan_beta_resummed(model),
               gm2calc::calculate_amu_2loop_non_tan_beta_resummed(mcpp));
}


TEST_CASE("uncertainty")
{
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   Cleanup_on_destruction cleanup(model);

   setup_gm2calc_scheme(model);

   const gm2calc::MSSMNoFV_onshell mcpp(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));

   CHECK_EQUAL(gm2calc_mssmnofv_calculate_uncertainty_amu_2loop(model),
               gm2calc::calculate_uncertainty_amu_2loop(mcpp));
}


TEST_CASE("uncertainty-size")
{
   MSSMNoFV_onshell* model = gm2calc_mssmnofv_new();
   Cleanup_on_destruction cleanup(model);

   setup_gm2calc_scheme(model);

   const gm2calc::MSSMNoFV_onshell mcpp(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));

   const auto damu_0l = gm2calc_mssmnofv_calculate_uncertainty_amu_0loop(model);
   const auto damu_1l = gm2calc_mssmnofv_calculate_uncertainty_amu_1loop(model);
   const auto damu_2l = gm2calc_mssmnofv_calculate_uncertainty_amu_2loop(model);

   CHECK_GT(damu_0l, damu_1l);
   CHECK_GT(damu_1l, damu_2l);
}
