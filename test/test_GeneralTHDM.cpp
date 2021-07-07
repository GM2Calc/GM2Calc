#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/GeneralTHDM.hpp"

#include <cmath>
#include <iostream>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


namespace {

double sqr(double x) noexcept { return x*x; }

const double pi = 3.1415926535897932;
const double alpha = 1./137.036;
const double e = std::sqrt(alpha*4*pi);
const double alpha_s = 0.1184;
const double g3 = std::sqrt(alpha_s*4*pi);
const double mw = 80.379;
const double mz = 91.1876;
const double cw = mw/mz;
const double sw = std::sqrt(1 - sqr(cw));
const double gY = e/cw;
const double g1 = gY*std::sqrt(5.0/3.0);
const double g2 = e/sw;
const double v = 2*mw/g2;

const Eigen::Matrix<double,3,3> ml{
   (Eigen::Matrix<double,3,3>()
    << 510.999e-6, 0.0, 0.0,
    0.0, 0.10565837, 0.0,
    0.0, 0.0, 1.7768).finished()
};

gm2calc::GeneralTHDM setup()
{
   const double tb = 20.0;
   const double sb = tb/std::sqrt(1.0 + sqr(tb));
   const double cb = 1./std::sqrt(1.0 + sqr(tb));
   const double vu = v*sb;
   const double vd = v*cb;

   gm2calc::GeneralTHDM model;
   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_v1(vd);
   model.set_v2(vu);
   model.set_Lambda1(0.7);
   model.set_Lambda2(0.6);
   model.set_Lambda3(0.5);
   model.set_Lambda4(0.4);
   model.set_Lambda5(0.3);
   model.set_Lambda6(0.2);
   model.set_Lambda7(0.1);
   model.set_M122(sqr(100.0));
   model.set_Ye(std::sqrt(2.0)*ml/vd);

   model.solve_ewsb();
   model.calculate_MSbar_masses();

   return model;
}

} // anonymous namespace


TEST_CASE("tree-level-spectrum")
{
   auto model = setup();
   const bool have_problem = model.get_problems().have_problem();

   CHECK(!have_problem);
}
