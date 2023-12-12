#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_linalg.hpp"
#include <limits>
#include <cmath>
#include <complex>


#define CHECK_SMALL(a,eps)                              \
   do {                                                 \
      CHECK(std::abs(a) < (eps));                       \
   } while (0)


template<class R_, class S_, int M_, int N_ = M_>
struct Test_fs {
    typedef R_ R;
    typedef S_ S;
    enum { M = M_ };
    enum { N = N_ };
};


TEST_CASE_TEMPLATE("test_casting_fs_svd", T,
    Test_fs<double, double, 1>,
    Test_fs<double, double, 6>,
    Test_fs<double, double, 4, 6>,
    Test_fs<double, double, 6, 4>,
    Test_fs<long double, long double, 1>,
    Test_fs<long double, long double, 6>,
    Test_fs<long double, long double, 4, 6>,
    Test_fs<long double, long double, 6, 4>
)
{
   typedef typename T::R R;
   constexpr Eigen::Index M = T::M;
   constexpr Eigen::Index N = T::N;
   const R eps = std::numeric_limits<R>::epsilon();

   Eigen::Matrix<R, M, N> m = Eigen::Matrix<R, M, N>::Random();
   Eigen::Array<R, MIN_(M, N), 1> s;
   Eigen::Matrix<std::complex<R>, M, M> u;
   Eigen::Matrix<std::complex<R>, N, N> v;

   gm2calc::fs_svd(m, s, u, v); // following SARAH convention
   Eigen::Matrix<std::complex<R>, M, N> sigma = u.conjugate() * m * v.adjoint();

   CHECK((s >= 0).all());
   for (Eigen::Index i = 0; i < sigma.rows(); i++) {
      for (Eigen::Index j = 0; j < sigma.cols(); j++) {
         CHECK_SMALL(std::abs(sigma(i, j) - (i == j ? s(i) : 0)), 50 * eps);
      }
   }

   for (Eigen::Index i = 0; i < s.size() - 1; i++) {
      CHECK(s[i] <= s[i + 1]);
   }

   gm2calc::fs_svd(m, s);
   CHECK((s >= 0).all());
   for (Eigen::Index i = 0; i < sigma.rows(); i++) {
      for (Eigen::Index j = 0; j < sigma.cols(); j++) {
         CHECK_SMALL(std::abs(sigma(i, j) - (i == j ? s(i) : 0)), 50 * eps);
      }
   }
}
