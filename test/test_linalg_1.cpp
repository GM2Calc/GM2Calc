#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_linalg.hpp"
#include <cmath>
#include <complex>


#define CHECK_SMALL(a,eps)                              \
   do {                                                 \
      CHECK(std::abs(a) < (eps));                       \
   } while (0)


template <class S_, int M_, int N_,
          void fxn_(const Eigen::Matrix<S_, M_, N_>&,
                    Eigen::Array<double, MIN_(M_, N_), 1>&,
                    Eigen::Matrix<S_, M_, M_>&, Eigen::Matrix<S_, N_, N_>&),
          void svs_(const Eigen::Matrix<S_, M_, N_>&,
                    Eigen::Array<double, MIN_(M_, N_), 1>&),
          bool check_ascending_order_ = false>
struct Test_svd {
   typedef S_ S;
   enum { M = M_ };
   enum { N = N_ };
   enum { check_ascending_order = check_ascending_order_ };
   void fxn(const Eigen::Matrix<S_, M_, N_>& m,
            Eigen::Array<double, MIN_(M_, N_), 1>& s,
            Eigen::Matrix<S_, M_, M_>& u, Eigen::Matrix<S_, N_, N_>& vh)
   {
      fxn_(m, s, u, vh);
   }
   void svs(const Eigen::Matrix<S_, M_, N_>& m,
            Eigen::Array<double, MIN_(M_, N_), 1>& s)
   {
      svs_(m, s);
   }
};


TEST_CASE_TEMPLATE("test_svd", T,
    Test_svd<std::complex<double>, 1, 1, gm2calc::svd, gm2calc::svd>,
    Test_svd<std::complex<double>, 2, 2, gm2calc::svd, gm2calc::svd>,
    Test_svd<std::complex<double>, 3, 3, gm2calc::svd, gm2calc::svd>,
    Test_svd<std::complex<double>, 6, 6, gm2calc::svd, gm2calc::svd>,
    Test_svd<double	         , 1, 1, gm2calc::svd, gm2calc::svd>,
    Test_svd<double	         , 2, 2, gm2calc::svd, gm2calc::svd>,
    Test_svd<double	         , 3, 3, gm2calc::svd, gm2calc::svd>,
    Test_svd<double	         , 6, 6, gm2calc::svd, gm2calc::svd>,
    Test_svd<std::complex<double>, 1, 1, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<std::complex<double>, 2, 2, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<std::complex<double>, 3, 3, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<std::complex<double>, 6, 6, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<std::complex<double>, 4, 6, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<std::complex<double>, 6, 4, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<double	         , 1, 1, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<double	         , 2, 2, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<double	         , 3, 3, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<double	         , 6, 6, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<double	         , 4, 6, gm2calc::reorder_svd, gm2calc::reorder_svd, true>,
    Test_svd<double	         , 6, 4, gm2calc::reorder_svd, gm2calc::reorder_svd, true>
   )
{
   typedef typename T::S S;
   constexpr Eigen::Index M = T::M;
   constexpr Eigen::Index N = T::N;

   const Eigen::Matrix<S, M, N> m = Eigen::Matrix<S, M, N>::Random();
   Eigen::Array<double, MIN_(M, N), 1> s;
   Eigen::Matrix<S, M, M> u;
   Eigen::Matrix<S, N, N> vh;

   T().fxn(m, s, u, vh);	// following LAPACK convention
   const Eigen::Matrix<S, M, N> sigma = u.adjoint() * m * vh.adjoint();

   CHECK((s >= 0).all());
   for (Eigen::Index i = 0; i < sigma.rows(); i++) {
      for (Eigen::Index j = 0; j < sigma.cols(); j++) {
         CHECK_SMALL(std::abs(sigma(i,j) - (i==j ? s(i) : 0)), 1e-13);
      }
   }

   if (T::check_ascending_order) {
      for (Eigen::Index i = 0; i < s.size()-1; i++) {
         CHECK(s[i] <= s[i+1]);
      }
   }

   T().svs(m, s);
   CHECK((s >= 0).all());
   for (Eigen::Index i = 0; i < sigma.rows(); i++) {
      for (Eigen::Index j = 0; j < sigma.cols(); j++) {
         CHECK_SMALL(std::abs(sigma(i,j) - (i==j ? s(i) : 0)), 1e-13);
      }
   }
}
