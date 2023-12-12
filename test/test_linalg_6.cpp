#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_linalg.hpp"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <limits>


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


template <typename R, typename S, Eigen::Index N>
void check_fs_diagonalize_symmetric(Eigen::Matrix<S, N, N> m)
{
   const R eps = std::numeric_limits<R>::epsilon();

   Eigen::Array<R, N, 1> s;
   Eigen::Matrix<std::complex<R>, N, N> u;

   gm2calc::fs_diagonalize_symmetric(m, s, u);
   const Eigen::Matrix<std::complex<R>, N, N> diag = u.conjugate() * m * u.adjoint();

   CHECK((s >= 0).all());
   for (Eigen::Index i = 0; i < N; i++) {
      for (Eigen::Index j = 0; j < N; j++) {
         CHECK_SMALL(std::abs(diag(i, j) - (i == j ? s(i) : 0)), 5000 * eps);
      }
   }

   for (Eigen::Index i = 0; i < N - 1; i++) {
      CHECK(s[i] <= s[i + 1]);
   }

   gm2calc::fs_diagonalize_symmetric(m, s);
   CHECK((s >= 0).all());
   for (Eigen::Index i = 0; i < N; i++) {
      for (Eigen::Index j = 0; j < N; j++) {
         CHECK_SMALL(std::abs(diag(i, j) - (i == j ? s(i) : 0)), 5000 * eps);
      }
   }
}


template <typename R>
R random_angle()
{
   return 2 * M_PI * std::rand() / RAND_MAX;
}


template<typename S, int N>
struct RandomUnitary;


template <typename R>
struct RandomUnitary<std::complex<R>, 1> {
   Eigen::Matrix<std::complex<R>, 1, 1> operator()() const
   {
      return Eigen::Matrix<std::complex<R>, 1, 1>(
         std::polar(R(1), random_angle<R>()));
   }
};


template <typename R>
struct RandomUnitary<R, 1> {
   Eigen::Matrix<R, 1, 1> operator()() const
   {
      return Eigen::Matrix<R, 1, 1>(2 * (rand() % 2) - 1);
   }
};


template<typename R>
struct RandomUnitary<std::complex<R>, 2> {
    Eigen::Matrix<std::complex<R>, 2, 2> operator()() const {
	R b = random_angle<R>();
	R c = random_angle<R>();
	R d = random_angle<R>();
	std::complex<R> p = std::polar(R(1), c);
	std::complex<R> q = std::polar(R(1), d);
	Eigen::Matrix<std::complex<R>, 2, 2> u;
	u <<           p *std::cos(b),           q *std::sin(b),
	    -std::conj(q)*std::sin(b), std::conj(p)*std::cos(b);
	u *= RandomUnitary<std::complex<R>, 1>()()(0,0);
	return u;
    }
};


template <typename R>
struct RandomUnitary<R, 2> {
   Eigen::Matrix<R, 2, 2> operator()() const
   {
      R b = random_angle<R>();
      Eigen::Matrix<R, 2, 2> o;
      o << std::cos(b), std::sin(b), -std::sin(b), std::cos(b);
      o *= RandomUnitary<R, 1>()()(0, 0);
      return o;
   }
};


template <typename S, int N>
struct RandomUnitary {
   Eigen::Matrix<S, N, N> operator()() const
   {
      Eigen::Matrix<S, N, N> u = Eigen::Matrix<S, N, N>::Identity();
      for (int i = 0; i < N - 1; i++) {
         for (int j = i + 1; j < N; j++) {
            Eigen::Matrix<S, N, N> s = Eigen::Matrix<S, N, N>::Identity();
            Eigen::Matrix<S, 2, 2> r = RandomUnitary<S, 2>()();
            s(i, i) = r(0, 0);
            s(i, j) = r(0, 1);
            s(j, i) = r(1, 0);
            s(j, j) = r(1, 1);
            u *= s;
         }
      }
      return u;
   }
};


TEST_CASE_TEMPLATE("test_fs_diagonalize_symmetric", T,
    // use Eigen::JacobiSVD
    Test_fs<double, std::complex<double>, 1>,
    Test_fs<double, std::complex<double>, 3>,
    Test_fs<double, std::complex<double>, 6>,
    Test_fs<long double, std::complex<long double>, 1>,
    Test_fs<long double, std::complex<long double>, 3>,
    Test_fs<long double, std::complex<long double>, 6>,
    // use Eigen::SelfAdjointEigenSolver
    Test_fs<double, double, 6>,
    Test_fs<long double, long double, 6>
)
{
   typedef typename T::R R;
   typedef typename T::S S;
   constexpr Eigen::Index N = T::N;

   Eigen::Matrix<S, N, N> m = Eigen::Matrix<S, N, N>::Random();
   m = ((m + m.transpose()) / 2).eval();
   check_fs_diagonalize_symmetric<R, S, N>(m);

   const Eigen::Matrix<S, N, N> u = RandomUnitary<S, N>()();
   const Eigen::Array<R, N, 1> s = Eigen::Array<R, N, 1>::Constant(2);
   m = u.transpose() * s.matrix().asDiagonal() * u;
   check_fs_diagonalize_symmetric<R, S, N>(m);
}
