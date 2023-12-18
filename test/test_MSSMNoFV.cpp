#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/MSSMNoFV_onshell.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2_linalg.hpp"
#include <iostream>
#include <Eigen/Core>


#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


namespace {

template <class Derived>
unsigned find_right_like_smuon(const Eigen::MatrixBase<Derived>& ZM)
{
   // Incorrect:
   // return (ZM(0,0) > ZM(0,1)) ? 1 : 0;

   return (std::abs(ZM(0,0)) > std::abs(ZM(0,1))) ? 1 : 0;
}

gm2calc::MSSMNoFV_onshell setup_slha(double eps)
{
   gm2calc::MSSMNoFV_onshell model;
   model.set_verbose_output(true);

   const double Pi = 3.141592653589793;
   const Eigen::Matrix<double,3,3> UnitMatrix
      = Eigen::Matrix<double,3,3>::Identity();

   // fill SM parameters
   model.set_alpha_MZ(0.0077552);               // 1L
   model.set_alpha_thompson(0.00729735);        // 2L
   model.set_g3(std::sqrt(4 * Pi * 0.1184));    // 2L
   model.get_physical().MFt   = 173.34;         // 2L
   model.get_physical().MFb   = 4.18;           // 2L, mb(mb) MS-bar
   model.get_physical().MFm   = 0.1056583715;   // 1L
   model.get_physical().MFtau = 1.777;          // 2L
   model.get_physical().MVWm  = 80.385;         // 1L
   model.get_physical().MVZ   = 91.1876;        // 1L

   // fill pole masses
   model.get_physical().MSvmL   =  5.18860573e+02; // 1L
   model.get_physical().MSm(0)  =  5.05095249e+02; // 1L
   model.get_physical().MSm(1)  =  5.25187016e+02; // 1L
   model.get_physical().MChi(0) =  2.01611468e+02; // 1L
   model.get_physical().MChi(1) =  4.10040273e+02; // 1L
   model.get_physical().MChi(2) =  5.16529941e+02; // 1L
   model.get_physical().MChi(3) =  5.45628749e+02; // 1L
   model.get_physical().MCha(0) =  4.09989890e+02; // 1L
   model.get_physical().MCha(1) =  5.46057190e+02; // 1L
   model.get_physical().MAh(1)  =  1.50000000e+03; // 2L

   // fill DR-bar parameters
   model.set_TB(40);                        // 1L
   model.set_Mu(500);                       // initial guess
   model.set_MassB(200);                    // initial guess
   model.set_MassWB(400);                   // initial guess
   model.set_MassG(2000);                   // 2L
   model.set_mq2(7000 * 7000 * UnitMatrix); // 2L
   model.set_ml2(0, 0, 500 * 500);          // 2L
   model.set_ml2(1, 1, 500 * 500);          // irrelevant
   model.set_ml2(2, 2, 500 * 500);          // 2L
   model.set_md2(7000 * 7000 * UnitMatrix); // 2L
   model.set_mu2(7000 * 7000 * UnitMatrix); // 2L
   model.set_me2(0, 0, 500 * 500);          // 2L
   model.set_me2(1, 1, 500 * 500);          // initial guess
   model.set_me2(2, 2, 500 * 500);          // 2L
   model.set_Au(2, 2, 0);                   // 2L
   model.set_Ad(2, 2, 0);                   // 2L
   model.set_Ae(1, 1, 0);                   // 1L
   model.set_Ae(2, 2, 0);                   // 2L
   model.set_scale(1000);                   // 2L

   // convert DR-bar parameters to on-shell
   model.convert_to_onshell(eps);

   // check for warnings
   if (model.get_problems().have_warning()) {
      std::cout << model.get_problems().get_warnings()  << '\n';
   }

   // check for problems
   if (model.get_problems().have_problem()) {
      std::cout << model.get_problems().get_problems()  << '\n';
   }

   return model;
}

gm2calc::MSSMNoFV_onshell setup_gm2calc()
{
   gm2calc::MSSMNoFV_onshell model;
   model.set_verbose_output(true);

   const double Pi = 3.141592653589793;
   const Eigen::Matrix<double,3,3> UnitMatrix
      = Eigen::Matrix<double,3,3>::Identity();

   // fill SM parameters
   model.set_alpha_MZ(0.0077552);               // 1L
   model.set_alpha_thompson(0.00729735);        // 2L
   model.set_g3(std::sqrt(4 * Pi * 0.1184));    // 2L
   model.get_physical().MFt   = 173.34;         // 2L
   model.get_physical().MFb   = 4.18;           // 2L, mb(mb) MS-bar
   model.get_physical().MFm   = 0.1056583715;   // 1L
   model.get_physical().MFtau = 1.777;          // 2L
   model.get_physical().MVWm  = 80.385;         // 1L
   model.get_physical().MVZ   = 91.1876;        // 1L

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

} // anonymous namespace


TEST_CASE("right_like_smuon")
{
   Eigen::Matrix<double,2,2> M;
   Eigen::Matrix<double,2,2> Z;
   Eigen::Array<double,2,1> m;

   M(0,0) = 1.0;
   M(0,1) = 0.1;
   M(1,0) = M(0,1);
   M(1,1) = 2.0;

   gm2calc::fs_diagonalize_hermitian<double,double,2>(M, m, Z);

   const auto idx = find_right_like_smuon(Z);
   const auto mR = m(idx); // right-like smuoon

   INFO("Z = " << Z);
   INFO("m = " << m.transpose());
   INFO("Mass of right-like smuon [" << idx << "] = " << mR);

   CHECK(std::abs(mR - M(1,1)) <= std::abs(mR - M(0,0)));
}


TEST_CASE("conversion_to_onshell")
{
   const double eps = 1e-8;

   gm2calc::MSSMNoFV_onshell model(setup_slha(eps));
   model.calculate_DRbar_masses();

   // Note: The parameters me2(1,1), ml2(1,1), Mu, MassB, MassWB are fixed by
   //
   // ml2(1,1) <-> MSvmL
   // me2(1,1) <-> MSm(0) (right-like smuon for this parameter point)
   // MassB    <-> MChi(0) (bino-like neutralino for this parameter point)
   // MassWB   <-> MCha(0) (wino-like chargino)
   // Mu       <-> MCha(1) (higgsino-like chargino)
   //
   // As a result, these masses will agree with the input pole masses
   // to the given precision.
   //
   // The other masses, MChi(1,2,3) and MSm(1), however, will not
   // necessarily agree to that precision, because they are
   // automatically fixed by all other parameters, in particular the
   // on-shell gauge couplings, the resummed Yukawa couplings.  In
   // other words, the other masses are no longer free parameters.
   // They can differ from the given input values because they are
   // calculated with other values for the gauge and Yukawa couplings.

   CHECK_CLOSE(model.get_MSvmL(), model.get_physical().MSvmL  , eps);
   CHECK_CLOSE(model.get_MSm(0) , model.get_physical().MSm(0) , eps);
   CHECK_CLOSE(model.get_MSm(1) , model.get_physical().MSm(1) , 1e-3);
   CHECK_CLOSE(model.get_MChi(0), model.get_physical().MChi(0), eps);
   CHECK_CLOSE(model.get_MChi(1), model.get_physical().MChi(1), 1e-3);
   CHECK_CLOSE(model.get_MChi(2), model.get_physical().MChi(2), 1e-2);
   CHECK_CLOSE(model.get_MChi(3), model.get_physical().MChi(3), 1e-3);
   CHECK_CLOSE(model.get_MCha(0), model.get_physical().MCha(0), eps);
   CHECK_CLOSE(model.get_MCha(1), model.get_physical().MCha(1), eps);
   CHECK_CLOSE(model.get_MA0()  , model.get_physical().MAh(1) , eps);

   CHECK_CLOSE(model.get_MFu()  , model.get_physical().MFu    , eps);
   CHECK_CLOSE(model.get_MFd()  , model.get_physical().MFd    , eps);
   CHECK_CLOSE(model.get_MFc()  , model.get_physical().MFc    , eps);
   CHECK_CLOSE(model.get_MFs()  , model.get_physical().MFs    , eps);
   // model.get_MFb() <-> mb(DR,Q), including higher-order corrections from the  tan(beta) resummation
   // physical.MFb    <-> mb(mb)
   // CHECK_CLOSE(model.get_MFb()  , model.get_physical().MFb    , eps);
   CHECK_CLOSE(model.get_MFt()  , model.get_physical().MFt    , eps);
   CHECK_CLOSE(model.get_MFe()  , model.get_physical().MFe    , eps);
   // model.get_MFm() includes higher-order corrections from the tan(beta) resummation
   CHECK_CLOSE(model.get_MFm()  , model.get_physical().MFm    , 0.01);
   // model.get_MFtau() includes higher-order corrections from the tan(beta) resummation
   CHECK_CLOSE(model.get_MFtau(), model.get_physical().MFtau  , 0.1);
}


TEST_CASE("calculate_masses")
{
   const double eps = 1e-15;

   gm2calc::MSSMNoFV_onshell model(setup_gm2calc());
   model.calculate_DRbar_masses();

   CHECK_CLOSE(model.get_MSvmL(), model.get_physical().MSvmL  , eps);
   CHECK_CLOSE(model.get_MSm(0) , model.get_physical().MSm(0) , eps);
   CHECK_CLOSE(model.get_MSm(1) , model.get_physical().MSm(1) , eps);
   CHECK_CLOSE(model.get_MChi(0), model.get_physical().MChi(0), eps);
   CHECK_CLOSE(model.get_MChi(1), model.get_physical().MChi(1), eps);
   CHECK_CLOSE(model.get_MChi(2), model.get_physical().MChi(2), eps);
   CHECK_CLOSE(model.get_MChi(3), model.get_physical().MChi(3), eps);
   CHECK_CLOSE(model.get_MCha(0), model.get_physical().MCha(0), eps);
   CHECK_CLOSE(model.get_MCha(1), model.get_physical().MCha(1), eps);
   CHECK_CLOSE(model.get_MA0()  , model.get_physical().MAh(1) , eps);

   CHECK_CLOSE(model.get_MFu()  , model.get_physical().MFu    , eps);
   CHECK_CLOSE(model.get_MFd()  , model.get_physical().MFd    , eps);
   CHECK_CLOSE(model.get_MFc()  , model.get_physical().MFc    , eps);
   CHECK_CLOSE(model.get_MFs()  , model.get_physical().MFs    , eps);
   // model.get_MFb() <-> mb(DR,Q), including higher-order corrections from the  tan(beta) resummation
   // physical.MFb    <-> mb(mb)
   // CHECK_CLOSE(model.get_MFb()  , model.get_physical().MFb    , eps);
   CHECK_CLOSE(model.get_MFt()  , model.get_physical().MFt    , eps);
   CHECK_CLOSE(model.get_MFe()  , model.get_physical().MFe    , eps);
   // model.get_MFm() includes higher-order corrections from the tan(beta) resummation
   CHECK_CLOSE(model.get_MFm()  , model.get_physical().MFm    , 0.01);
   // model.get_MFtau() includes higher-order corrections from the tan(beta) resummation
   CHECK_CLOSE(model.get_MFtau(), model.get_physical().MFtau  , 0.1);
}


// This test ensures that the Goldstone bosons are always at the first
// position (index 0) in the MAh and MHpm multiplets, even if they are
// heavier than the corresponding Higgs bosons.
TEST_CASE("goldstone_boson_position")
{
   const double eps = 1e-15;
   const double Pi = 3.141592653589793;
   const Eigen::Matrix<double,3,3> UnitMatrix
      = Eigen::Matrix<double,3,3>::Identity();

   gm2calc::MSSMNoFV_onshell model;

   model.set_alpha_MZ(0.0077552);
   model.set_alpha_thompson(0.00729735);
   model.set_g3(std::sqrt(4 * Pi * 0.1184));
   model.get_physical().MFt   = 173.34;
   model.get_physical().MFb   = 4.18;
   model.get_physical().MFm   = 0.1056583715;
   model.get_physical().MFtau = 1.777;
   model.get_physical().MVWm  = 80.385;
   model.get_physical().MVZ   = 91.1876;

   const double cW = model.get_physical().MVWm / model.get_physical().MVZ;
   model.set_g1(std::sqrt(5./3.) * model.get_EL() / cW);
   model.set_g2(model.get_EL()/std::sqrt(1. - cW*cW));
   model.set_TB(2);
   model.set_Mu(500);
   model.set_MassB(200);
   model.set_MassWB(400);
   model.set_MassG(2000);
   model.set_BMu(500); // very small so that MAh < MZ
   model.set_Mu(350);
   model.set_MassB(150);
   model.set_MassWB(300);
   model.set_MassG(1000);
   model.set_mq2(500 * 500 * UnitMatrix);
   model.set_ml2(500 * 500 * UnitMatrix);
   model.set_md2(500 * 500 * UnitMatrix);
   model.set_mu2(500 * 500 * UnitMatrix);
   model.set_me2(500 * 500 * UnitMatrix);
   model.set_Au(2,2,0);
   model.set_Ad(2,2,0);
   model.set_Ae(2,2,0);
   model.set_scale(1000);

   model.calculate_DRbar_masses();

   CHECK(!model.get_problems().have_problem());

   CHECK_CLOSE(model.get_MAh(0), model.get_MVZ(), eps);
   CHECK_CLOSE(model.get_MHpm(0), model.get_MVWm(), eps);
}


TEST_CASE("uncertainty-size")
{
   const auto model = setup_gm2calc();
   const auto damu_0l = gm2calc::calculate_uncertainty_amu_0loop(model);
   const auto damu_1l = gm2calc::calculate_uncertainty_amu_1loop(model);
   const auto damu_2l = gm2calc::calculate_uncertainty_amu_2loop(model);

   CHECK_GT(damu_0l, damu_1l);
   CHECK_GT(damu_1l, damu_2l);
}


TEST_CASE("print")
{
   const auto model = setup_gm2calc();
   std::cout << model;
   std::cout << model.get_physical();
}
