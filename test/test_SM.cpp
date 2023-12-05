#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/SM.hpp"
#include <cmath>
#include <limits>


#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


bool is_unitary(const Eigen::Matrix<std::complex<double>,3,3>& m, double eps)
{
   const Eigen::Matrix<std::complex<double>,3,3> unit = Eigen::Matrix<std::complex<double>,3,3>::Identity();
   const Eigen::Matrix<std::complex<double>,3,3> cca = m * m.adjoint();
   const double max_diff = (cca - unit).cwiseAbs().maxCoeff();
   return max_diff <= eps;
}


// check unitarity of default CKM matrix
TEST_CASE("test-CKM-unitarity")
{
   const double eps = std::numeric_limits<double>::epsilon();
   CHECK(is_unitary(gm2calc::SM().get_ckm(), eps));
}


// check unitarity of CKM matrix constructed from angles
TEST_CASE("test-CKM-unitarity")
{
   const double eps = std::numeric_limits<double>::epsilon();

   auto sm = gm2calc::SM();
   sm.set_ckm_from_angles(0.1, 0.2, 0.3, 0.4);

   CHECK(is_unitary(sm.get_ckm(), eps));
}


// calculation of e_0
TEST_CASE("test-e0")
{
   const double eps = std::numeric_limits<double>::epsilon();
   const double pi = 3.1415926535897932;
   const double e0 = 1.0/127.0;
   const double alpha_em_0 = e0*e0/(4*pi);

   auto sm = gm2calc::SM();
   sm.set_alpha_em_0(alpha_em_0);

   CHECK_CLOSE(sm.get_e_0(), e0, eps);
}


// calculation of e(MZ)
TEST_CASE("test-eMZ")
{
   const double eps = std::numeric_limits<double>::epsilon();
   const double pi = 3.1415926535897932;
   const double eMZ = 1.0/127.0;
   const double alpha_em_mz = eMZ*eMZ/(4*pi);

   auto sm = gm2calc::SM();
   sm.set_alpha_em_mz(alpha_em_mz);

   CHECK_CLOSE(sm.get_e_mz(), eMZ, eps);
}


// calculation of g3
TEST_CASE("test-g3")
{
   const double eps = std::numeric_limits<double>::epsilon();
   const double pi = 3.1415926535897932;
   const double g3 = 0.1;
   const double alpha_s_mz = g3*g3/(4*pi);

   auto sm = gm2calc::SM();
   sm.set_alpha_s_mz(alpha_s_mz);

   CHECK_CLOSE(sm.get_g3(), g3, eps);
}


// calculation of gY and g2
TEST_CASE("test-g3")
{
   const double eps = std::numeric_limits<double>::epsilon();
   const double pi = 3.1415926535897932;
   const auto sm = gm2calc::SM();
   const double gY = sm.get_gY();
   const double g2 = sm.get_g2();
   const double eMZ = sm.get_e_mz();
   const double theta = std::acos(sm.get_mw()/sm.get_mz());

   CHECK_CLOSE(std::cos(theta), sm.get_cw(), eps);
   CHECK_CLOSE(std::sin(theta), sm.get_sw(), eps);
   CHECK_CLOSE(eMZ, gY*std::cos(theta), eps);
   CHECK_CLOSE(eMZ, g2*std::sin(theta), eps);
}
