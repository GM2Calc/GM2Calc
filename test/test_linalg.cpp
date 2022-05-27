#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_linalg.hpp"
#include <complex>
#include <limits>
#include <Eigen/Core>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


using CM22 = Eigen::Matrix<std::complex<double>,2,2>;
using RM22 = Eigen::Matrix<double,2,2>;
using RA21 = Eigen::Array<double,2,1>;


RM22 create_real_symmetric_2x2(double a, double b, double c)
{
   return (RM22() << a, b, b, c).finished();
}


// diagonalization using the general SelfAdjointEigenSolver
template<class Real, class Scalar, int>
void fs_diagonalize_hermitian_2x2
(const Eigen::Matrix<Scalar, 2, 2>& m,
 Eigen::Array<Real, 2, 1>& w,
 Eigen::Matrix<Scalar, 2, 2>& z)
{
   constexpr int N = 2;

   Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar,N,N> >
      es(m, Eigen::ComputeEigenvectors);
   w = es.eigenvalues();
   z = es.eigenvectors();

   Eigen::PermutationMatrix<N> p;
   p.setIdentity();
   std::sort(p.indices().data(), p.indices().data() + p.indices().size(),
             [&w] (int i, int j) { return std::abs(w[i]) < std::abs(w[j]); });
#if EIGEN_VERSION_AT_LEAST(3,1,4)
   w.matrix().transpose() *= p;
#else
   Eigen::Map<Eigen::Matrix<Real, N, 1> >(w.data()).transpose() *= p;
#endif
   z = (z * p).adjoint().eval();
}


template<class Real, class Scalar>
void test_hermitian_2x2(const Eigen::Matrix<Scalar, 2, 2>& m, double eps)
{
   Eigen::Array<Real, 2, 1> v1, v2;
   Eigen::Matrix<Scalar, 2, 2> z1, z2;

   fs_diagonalize_hermitian_2x2<Real, Scalar, 2>(m, v1, z1);
   gm2calc::fs_diagonalize_hermitian<Real, Scalar, 2>(m, v2, z2);

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

   struct Matrix_elements {
      double a, b, c, eps;
   } const data[] = {
      { 1, 0,  2, 10*eps},
      { 2, 0,  1, 10*eps},
      { 2, 5,  1, 10*eps},
      {-1, 0, -2, 10*eps},
      {-2, 0, -1, 10*eps},
      {-2, 5, -1, 10*eps},
      { 1, 2,  3, 10*eps}
   };

   for (const auto& d: data) {
      test_hermitian_2x2<double, double>(create_real_symmetric_2x2(d.a, d.b, d.c), d.eps);
   }
}
