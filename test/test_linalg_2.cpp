#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_linalg.hpp"
#include <cmath>
#include <complex>


#define CHECK_SMALL(a,eps)                              \
   do {                                                 \
      CHECK(std::abs(a) < (eps));                       \
   } while (0)


template <
   class S_, int N_,
   void fxn_(const Eigen::Matrix<S_, N_, N_>&, Eigen::Array<double, N_, 1>&,
             Eigen::Matrix<std::complex<double>, N_, N_>&),
   void svs_(const Eigen::Matrix<S_, N_, N_>&, Eigen::Array<double, N_, 1>&),
   bool check_ascending_order_ = false>
struct Test_diagonalize_symmetric {
   typedef S_ S;
   enum { N = N_ };
   enum { check_ascending_order = check_ascending_order_ };
   void fxn(const Eigen::Matrix<S_, N_, N_>& m, Eigen::Array<double, N_, 1>& s,
            Eigen::Matrix<std::complex<double>, N_, N_>& u)
   {
      fxn_(m, s, u);
   }
   void svs(const Eigen::Matrix<S_, N_, N_>& m, Eigen::Array<double, N_, 1>& s)
   {
      svs_(m, s);
   }
};

TEST_CASE_TEMPLATE("test_diagonalize_symmetric", T,
    // use Eigen::JacobiSVD
    Test_diagonalize_symmetric<std::complex<double>, 1, gm2calc::diagonalize_symmetric        , gm2calc::diagonalize_symmetric>,
    Test_diagonalize_symmetric<std::complex<double>, 2, gm2calc::diagonalize_symmetric        , gm2calc::diagonalize_symmetric>,
    Test_diagonalize_symmetric<std::complex<double>, 3, gm2calc::diagonalize_symmetric        , gm2calc::diagonalize_symmetric>,
    Test_diagonalize_symmetric<std::complex<double>, 4, gm2calc::diagonalize_symmetric        , gm2calc::diagonalize_symmetric>,
    Test_diagonalize_symmetric<std::complex<double>, 6, gm2calc::diagonalize_symmetric        , gm2calc::diagonalize_symmetric>,
    Test_diagonalize_symmetric<std::complex<double>, 1, gm2calc::reorder_diagonalize_symmetric, gm2calc::reorder_diagonalize_symmetric, true>,
    Test_diagonalize_symmetric<std::complex<double>, 2, gm2calc::reorder_diagonalize_symmetric, gm2calc::reorder_diagonalize_symmetric, true>,
    Test_diagonalize_symmetric<std::complex<double>, 3, gm2calc::reorder_diagonalize_symmetric, gm2calc::reorder_diagonalize_symmetric, true>,
    Test_diagonalize_symmetric<std::complex<double>, 4, gm2calc::reorder_diagonalize_symmetric, gm2calc::reorder_diagonalize_symmetric, true>,
    Test_diagonalize_symmetric<std::complex<double>, 6, gm2calc::reorder_diagonalize_symmetric, gm2calc::reorder_diagonalize_symmetric, true>,
    // use Eigen::SelfAdjointEigenSolver
    Test_diagonalize_symmetric<double, 6, gm2calc::diagonalize_symmetric        , gm2calc::diagonalize_symmetric>,
    Test_diagonalize_symmetric<double, 6, gm2calc::reorder_diagonalize_symmetric, gm2calc::reorder_diagonalize_symmetric, true>
)
{
   typedef typename T::S S;
   constexpr Eigen::Index N = T::N;

   Eigen::Matrix<S, N, N> m = Eigen::Matrix<S, N, N>::Random();
   m = ((m + m.transpose())/2).eval();
   Eigen::Array<double, N, 1> s;
   Eigen::Matrix<std::complex<double>, N, N> u;

   T().fxn(m, s, u);
   const Eigen::Matrix<std::complex<double>, N, N> diag = u.adjoint() * m * u.conjugate();

   CHECK((s >= 0).all());
   for (Eigen::Index i = 0; i < N; i++) {
      for (Eigen::Index j = 0; j < N; j++) {
         CHECK_SMALL(std::abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-12);
      }
   }

   if (T::check_ascending_order) {
      for (Eigen::Index i = 0; i < N-1; i++) {
         CHECK(s[i] <= s[i+1]);
      }
   }

   T().svs(m, s);
   CHECK((s >= 0).all());
   for (Eigen::Index i = 0; i < N; i++) {
      for (Eigen::Index j = 0; j < N; j++) {
         CHECK_SMALL(std::abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-12);
      }
   }
}
