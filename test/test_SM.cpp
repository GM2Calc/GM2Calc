#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/SM.hpp"
#include <limits>


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
