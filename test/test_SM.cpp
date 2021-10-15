#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/SM.hpp"
#include <limits>

// check unitarity of default CKM matrix
TEST_CASE("test-CKM-unitarity")
{
   const double eps = std::numeric_limits<double>::epsilon();

   const Eigen::Matrix<std::complex<double>,3,3> ckm = gm2calc::SM().get_ckm();
   const Eigen::Matrix<std::complex<double>,3,3> unit = Eigen::Matrix<std::complex<double>,3,3>::Identity();
   const Eigen::Matrix<std::complex<double>,3,3> cca = ckm * ckm.adjoint();
   const double max_diff = (cca - unit).cwiseAbs().maxCoeff();

   INFO("max diff: " << max_diff);
   INFO("eps: " << eps);

   CHECK(max_diff <= eps);
}
