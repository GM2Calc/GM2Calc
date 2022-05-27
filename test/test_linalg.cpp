#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_linalg.hpp"
#include <limits>
#include <Eigen/Core>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)

using RM22 = Eigen::Matrix<double,2,2>;
using AR21 = Eigen::Array<double,2,1>;

RM22 create_real_symmetric_2x2(double a, double b, double c)
{
   return (RM22() << a, b, b, c).finished();
}

void test_real_symmetric_2x2(const RM22& m, double eps)
{
   AR21 v1, v2;
   RM22 z1, z2;

   gm2calc::fs_diagonalize_hermitian<double,double,2>(m, v1, z1);
   gm2calc::fs_diagonalize_hermitian_2x2<double>(m, v2, z2);

   CHECK_CLOSE(v1(0), v2(0), eps);
   CHECK_CLOSE(v1(1), v2(1), eps);

   // check sorted
   CHECK(std::abs(v1(0)) <= std::abs(v1(1)));
   CHECK(std::abs(v2(0)) <= std::abs(v2(1)));

   // check that z diagonalizes m and yields v
   CHECK_CLOSE((m - z1.transpose() * v1.matrix().asDiagonal() * z1).cwiseAbs().maxCoeff(), 0.0, eps);
   CHECK_CLOSE((m - z2.transpose() * v2.matrix().asDiagonal() * z2).cwiseAbs().maxCoeff(), 0.0, eps);
}

// test eigenvalues of real symmetric 2x2 matrix
TEST_CASE("test-real-symmetric-2x2")
{
   const double eps = std::numeric_limits<double>::epsilon();

   test_real_symmetric_2x2(create_real_symmetric_2x2(1, 0, 2), 10*eps);
   test_real_symmetric_2x2(create_real_symmetric_2x2(2, 0, 1), 10*eps);
   test_real_symmetric_2x2(create_real_symmetric_2x2(2, 5, 1), 10*eps);
   test_real_symmetric_2x2(create_real_symmetric_2x2(-1, 0, -2), 10*eps);
   test_real_symmetric_2x2(create_real_symmetric_2x2(-2, 0, -1), 10*eps);
   test_real_symmetric_2x2(create_real_symmetric_2x2(-2, 5, -1), 10*eps);
   test_real_symmetric_2x2(create_real_symmetric_2x2(1, 2, 3), 10*eps);
}
