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

void test_eigensystem_real_symmetric_2x2(const RM22& m, double eps)
{
   AR21 v1, v2;
   RM22 z1, z2;

   gm2calc::fs_diagonalize_hermitian(m, v1, z1);
   gm2calc::fs_diagonalize_hermitian_2x2(m, v2, z2);

   CHECK_CLOSE(v1(0), v2(0), eps);
   CHECK_CLOSE(v1(1), v2(1), eps);
   CHECK_CLOSE(z1(0,0), z2(0,0), eps);
   CHECK_CLOSE(z1(0,1), z2(0,1), eps);
   CHECK_CLOSE(z1(1,0), z2(1,0), eps);
   CHECK_CLOSE(z1(1,1), z2(1,1), eps);
}

// test eigenvalues of real symmetric 2x2 matrix
TEST_CASE("test-real-symmetric-2x2")
{
   const double eps = std::numeric_limits<double>::epsilon();

   test_eigensystem_real_symmetric_2x2(create_real_symmetric_2x2(1.0, 2.0, 3.0), 10*eps);
}
