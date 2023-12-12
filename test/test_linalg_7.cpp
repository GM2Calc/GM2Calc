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


template <class R_, class S_, int M_, int N_ = M_>
struct Test_fs {
   typedef R_ R;
   typedef S_ S;
   enum { M = M_ };
   enum { N = N_ };
};


#ifdef TEST_LINALG2_PART7
using namespace boost::mpl::placeholders;

typedef boost::mpl::fold<
    boost::mpl::range_c<int, 0, 10>,
    boost::mpl::list<>,
    boost::mpl::push_front<
      boost::mpl::push_front<
        boost::mpl::push_front<
          boost::mpl::push_front<
	      _1,
	      boost::mpl::pair<Test_fs<double, std::complex<double>, 6>, _2> >,
	    boost::mpl::pair<Test_fs<double, double, 6>, _2> >,
	  boost::mpl::pair<Test_fs<long double, std::complex<long double>, 6>,_2> >,
	boost::mpl::pair<Test_fs<long double, long double, 6>, _2> >
>::type fs_diagonalize_hermitian_tests;


TEST_CASE_TEMPLATE("test_fs_diagonalize_hermitian", P,
                   fs_diagonalize_hermitian_tests)
{
   typedef typename P::first T;
   typedef typename T::R R;
   typedef typename T::S S;
   constexpr Eigen::Index N = T::N;
   const R eps = std::numeric_limits<R>::epsilon();

   Eigen::Matrix<S, N, N> m = Eigen::Matrix<S, N, N>::Random();
   m = ((m + m.adjoint()) / 2).eval();
   Eigen::Array<R, N, 1> w;
   Eigen::Matrix<S, N, N> z;

   gm2calc::fs_diagonalize_hermitian(m, w, z); // following SARAH convention
   const Eigen::Matrix<S, N, N> diag = z * m * z.adjoint();

   for (Eigen::Index i = 0; i < N; i++) {
      for (Eigen::Index j = 0; j < N; j++) {
         CHECK_SMALL(std::abs(diag(i, j) - (i == j ? w(i) : 0)), 50000 * eps);
      }
   }

   for (Eigen::Index i = 0; i < N - 1; i++) {
      CHECK(std::abs(w[i]) <= std::abs(w[i + 1]));
   }

   gm2calc::fs_diagonalize_hermitian(m, w);
   for (Eigen::Index i = 0; i < N; i++) {
      for (Eigen::Index j = 0; j < N; j++) {
         CHECK_SMALL(std::abs(diag(i, j) - (i == j ? w(i) : 0)), 50000 * eps);
      }
   }
}
#endif // TEST_LINALG2_PART7
