#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"

#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2calc/MSSMNoFV_onshell.hpp"

#include "stopwatch.hpp"

#include <random>
#include <utility>

namespace {

double sqr(double x) noexcept { return x*x; }

/// random number generator
double random(double start, double stop)
{
   static std::minstd_rand gen;
   std::uniform_real_distribution<double> dist(start, stop);
   return dist(gen);
}

/// random value for SUSY scale MS
auto rMS = [] {
   const double MS_min = 400;
   const double MS_max = 4000;
   return random(MS_min, MS_max);
};

/// random value for tan(beta)
auto rTB = [] { return random(2, 1000); };

/// generate random SLHA input parameter point
gm2calc::MSSMNoFV_onshell random_point_slha()
{
   gm2calc::MSSMNoFV_onshell model;

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
   model.get_physical().MSvmL   = rMS();    // 1L
   model.get_physical().MSm(0)  = rMS();    // 1L
   model.get_physical().MSm(1)  = rMS();    // 1L
   model.get_physical().MChi(0) = rMS();    // 1L
   model.get_physical().MChi(1) = rMS();    // 1L
   model.get_physical().MChi(2) = rMS();    // 1L
   model.get_physical().MChi(3) = rMS();    // 1L
   model.get_physical().MCha(0) = rMS();    // 1L
   model.get_physical().MCha(1) = rMS();    // 1L
   model.get_physical().MAh(1)  = rMS();    // 2L

   // fill DR-bar parameters
   model.set_TB(rTB());         // 1L
   model.set_Mu(rMS());                     // initial guess
   model.set_MassB(rMS());                  // initial guess
   model.set_MassWB(rMS());                 // initial guess
   model.set_MassG(rMS());                  // 2L
   model.set_mq2(sqr(rMS()) * UnitMatrix);  // 2L
   model.set_ml2(sqr(rMS()) * UnitMatrix);  // 2L
   model.set_md2(sqr(rMS()) * UnitMatrix);  // 2L
   model.set_mu2(sqr(rMS()) * UnitMatrix);  // 2L
   model.set_me2(sqr(rMS()) * UnitMatrix);  // 2L
   model.set_Au(0.0 * UnitMatrix);          // 2L
   model.set_Ad(0.0 * UnitMatrix);          // 2L
   model.set_Ae(0.0 * UnitMatrix);          // 2L
   model.set_Ae(1, 1, 0);                   // 1L
   model.set_scale(rMS());                  // 2L

   // convert DR-bar parameters to on-shell
   model.convert_to_onshell();

   return model;
}

/// generate random GM2Calc input parameter point
gm2calc::MSSMNoFV_onshell random_point_gm2calc()
{
   gm2calc::MSSMNoFV_onshell model;

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

   // fill SUSY parameters
   model.set_TB(rTB());                     // 1L
   model.set_Mu(rMS());                     // 1L
   model.set_MassB(rMS());                  // 1L
   model.set_MassWB(rMS());                 // 1L
   model.set_MassG(rMS());                  // 2L
   model.set_mq2(sqr(rMS()) * UnitMatrix);  // 2L
   model.set_ml2(sqr(rMS()) * UnitMatrix);  // 2L
   model.set_md2(sqr(rMS()) * UnitMatrix);  // 2L
   model.set_mu2(sqr(rMS()) * UnitMatrix);  // 2L
   model.set_me2(sqr(rMS()) * UnitMatrix);  // 2L
   model.set_Au(0.0 * UnitMatrix);          // 2L
   model.set_Ad(0.0 * UnitMatrix);          // 2L
   model.set_Ae(0.0 * UnitMatrix);          // 2L
   model.set_Ae(1, 1, 0);                   // 1L
   model.set_scale(rMS());                  // 2L
   model.set_MA0(rMS());                    // 2L

   // calculate mass spectrum
   model.calculate_masses();

   return model;
}

/// calculate amu and uncertainty
std::pair<double, double> calculate_amu(const gm2calc::MSSMNoFV_onshell& model)
{
   const double amu = gm2calc::calculate_amu_1loop(model) +
                      gm2calc::calculate_amu_2loop(model);

   const double delta_amu = gm2calc::calculate_uncertainty_amu_2loop(model);

   return std::make_pair(amu, delta_amu);
}

/// measure time for calling the function f()
template <class F>
double time_in_milliseconds(F&& f)
{
   gm2calc::Stopwatch sw;
   sw.start();
   f();
   sw.stop();
   return sw.get_time_in_milliseconds();
}

/// measure time for calling the function f() N times
template <class F>
double time_in_milliseconds(unsigned N, F&& f)
{
   auto loop = [N, f]() {
      for (unsigned i = 0; i < N; i++) {
         try {
            (void) f();
         } catch (...) {
         }
      }
   };

   return time_in_milliseconds(loop);
}

} // anonymous namespace

TEST_CASE("benchmark GM2Calc")
{
   const unsigned N = 10000;
   const auto time_in_ms = time_in_milliseconds(
      N, [] { return calculate_amu(random_point_gm2calc()); });

   std::cout << "Average time per point: " << time_in_ms/N << " ms (GM2Calc)\n";
}

TEST_CASE("benchmark SLHA")
{
   const unsigned N = 1000;
   const auto time_in_ms = time_in_milliseconds(
      N, [] { return calculate_amu(random_point_slha()); });

   std::cout << "Average time per point: " << time_in_ms/N << " ms (SLHA)\n";
}
