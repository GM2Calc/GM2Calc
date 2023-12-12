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


template <class S_, int N_,
          void fxn_(const Eigen::Matrix<S_, N_, N_>&,
                    Eigen::Array<double, N_, 1>&, Eigen::Matrix<S_, N_, N_>*)>
struct Test_diagonalize_hermitian {
   typedef S_ S;
   enum { N = N_ };
   void fxn(const Eigen::Matrix<S_, N_, N_>& m, Eigen::Array<double, N_, 1>& w,
            Eigen::Matrix<S_, N_, N_>* z)
   {
      fxn_(m, w, z);
   }
};


TEST_CASE_TEMPLATE("test_diagonalize_hermitian", T,
   // use Eigen::SelfAdjointEigenSolver
   Test_diagonalize_hermitian<std::complex<double>, 6, gm2calc::hermitian_eigen>,
   Test_diagonalize_hermitian<double              , 6, gm2calc::hermitian_eigen>
)
{
   typedef typename T::S S;
   constexpr Eigen::Index N = T::N;

   Eigen::Matrix<S, N, N> m = Eigen::Matrix<S, N, N>::Random();
   m = ((m + m.adjoint()) / 2).eval();
   Eigen::Array<double, N, 1> w;
   Eigen::Matrix<S, N, N> z;

   T().fxn(m, w, &z); // following LAPACK convention
   const Eigen::Matrix<S, N, N> diag = z.adjoint() * m * z;

   for (Eigen::Index i = 0; i < N; i++) {
      for (Eigen::Index j = 0; j < N; j++) {
         CHECK_SMALL(abs(diag(i, j) - (i == j ? w(i) : 0)), 1e-12);
      }
   }

   for (Eigen::Index i = 0; i < N - 1; i++) {
      CHECK(w[i] <= w[i + 1]);
   }
}
