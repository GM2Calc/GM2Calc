#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_linalg.hpp"
#include <limits>
#include <cmath>
#include <complex>
#include <iostream>

// #define BOOST_TEST_DYN_LINK
// #define BOOST_TEST_MODULE test_linalg2

// #include <boost/test/unit_test.hpp>
// #include <boost/mpl/pair.hpp>
// #include <boost/mpl/list.hpp>
// #include <boost/mpl/int.hpp>
// #include <boost/mpl/range_c.hpp>
// #include <boost/mpl/fold.hpp>
// #include <boost/mpl/push_front.hpp>
// #include <boost/mpl/placeholders.hpp>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)

#define CHECK_SMALL(a,eps)                              \
   do {                                                 \
      CHECK(std::abs(a) < (eps));                       \
   } while (0)


template<class S_, int M_, int N_,
	 void fxn_(const Eigen::Matrix<S_, M_, N_>&,
		   Eigen::Array<double, MIN_(M_, N_), 1>&,
		   Eigen::Matrix<S_, M_, M_>&,
		   Eigen::Matrix<S_, N_, N_>&),
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
	     Eigen::Matrix<S_, M_, M_>& u,
	     Eigen::Matrix<S_, N_, N_>& vh)
    { fxn_(m, s, u, vh); }
    void svs(const Eigen::Matrix<S_, M_, N_>& m,
	     Eigen::Array<double, MIN_(M_, N_), 1>& s)
    { svs_(m, s); }
};

#ifdef TEST_LINALG2_PART1
typedef boost::mpl::list<
    // use Eigen::JacobiSVD
    Test_svd<std::complex<double>, 1, 1, svd, svd>,
    Test_svd<std::complex<double>, 2, 2, svd, svd>,
    Test_svd<std::complex<double>, 3, 3, svd, svd>,
    Test_svd<std::complex<double>, 6, 6, svd, svd>,
    Test_svd<double	    , 1, 1, svd, svd>,
    Test_svd<double	    , 2, 2, svd, svd>,
    Test_svd<double	    , 3, 3, svd, svd>,
    Test_svd<double	    , 6, 6, svd, svd>,

    Test_svd<std::complex<double>, 1, 1, reorder_svd, reorder_svd, true>,
    Test_svd<std::complex<double>, 2, 2, reorder_svd, reorder_svd, true>,
    Test_svd<std::complex<double>, 3, 3, reorder_svd, reorder_svd, true>,
    Test_svd<std::complex<double>, 6, 6, reorder_svd, reorder_svd, true>,
    Test_svd<std::complex<double>, 4, 6, reorder_svd, reorder_svd, true>,
    Test_svd<std::complex<double>, 6, 4, reorder_svd, reorder_svd, true>,
    Test_svd<double	    , 1, 1, reorder_svd, reorder_svd, true>,
    Test_svd<double	    , 2, 2, reorder_svd, reorder_svd, true>,
    Test_svd<double	    , 3, 3, reorder_svd, reorder_svd, true>,
    Test_svd<double	    , 6, 6, reorder_svd, reorder_svd, true>,
    Test_svd<double	    , 4, 6, reorder_svd, reorder_svd, true>,
    Test_svd<double	    , 6, 4, reorder_svd, reorder_svd, true>
> svd_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_svd, T, svd_tests)
{
    typedef typename T::S S;
    const std::size_t M = T::M;
    const std::size_t N = T::N;

    Eigen::Matrix<S, M, N> m = Eigen::Matrix<S, M, N>::Random();
    Eigen::Array<double, MIN_(M, N), 1> s;
    Eigen::Matrix<S, M, M> u;
    Eigen::Matrix<S, N, N> vh;

    T().fxn(m, s, u, vh);	// following LAPACK convention
    Eigen::Matrix<S, M, N> sigma = u.adjoint() * m * vh.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (Eigen::Index i = 0; i < sigma.rows(); i++)
	for (Eigen::Index j = 0; j < sigma.cols(); j++)
	    BOOST_CHECK_SMALL(abs(sigma(i,j) - (i==j ? s(i) : 0)), 1e-13);

    if (T::check_ascending_order)
	for (Eigen::Index i = 0; i < s.size()-1; i++)
	    BOOST_CHECK(s[i] <= s[i+1]);

    T().svs(m, s);
    BOOST_CHECK((s >= 0).all());
    for (Eigen::Index i = 0; i < sigma.rows(); i++)
	for (Eigen::Index j = 0; j < sigma.cols(); j++)
	    BOOST_CHECK_SMALL(abs(sigma(i,j) - (i==j ? s(i) : 0)), 1e-13);
}
#endif // TEST_LINALG2_PART1

template<class S_, int N_,
	 void fxn_(const Eigen::Matrix<S_, N_, N_>&,
		   Eigen::Array<double, N_, 1>&, Eigen::Matrix<std::complex<double>, N_, N_>&),
	 void svs_(const Eigen::Matrix<S_, N_, N_>&,
		   Eigen::Array<double, N_, 1>&),
	 bool check_ascending_order_ = false>
struct Test_diagonalize_symmetric {
    typedef S_ S;
    enum { N = N_ };
    enum { check_ascending_order = check_ascending_order_ };
    void fxn(const Eigen::Matrix<S_, N_, N_>& m,
	     Eigen::Array<double, N_, 1>& s, Eigen::Matrix<std::complex<double>, N_, N_>& u)
    { fxn_(m, s, u); }
    void svs(const Eigen::Matrix<S_, N_, N_>& m,
	     Eigen::Array<double, N_, 1>& s)
    { svs_(m, s); }
};

#ifdef TEST_LINALG2_PART2
typedef boost::mpl::list<
    // use Eigen::JacobiSVD
    Test_diagonalize_symmetric
	<std::complex<double>, 1, diagonalize_symmetric, diagonalize_symmetric>,
    Test_diagonalize_symmetric
	<std::complex<double>, 2, diagonalize_symmetric, diagonalize_symmetric>,
    Test_diagonalize_symmetric
	<std::complex<double>, 3, diagonalize_symmetric, diagonalize_symmetric>,
    Test_diagonalize_symmetric
	<std::complex<double>, 4, diagonalize_symmetric, diagonalize_symmetric>,
    Test_diagonalize_symmetric
	<std::complex<double>, 6, diagonalize_symmetric, diagonalize_symmetric>,

    Test_diagonalize_symmetric
	<std::complex<double>, 1,
	 reorder_diagonalize_symmetric, reorder_diagonalize_symmetric, true>,
    Test_diagonalize_symmetric
	<std::complex<double>, 2,
	 reorder_diagonalize_symmetric, reorder_diagonalize_symmetric, true>,
    Test_diagonalize_symmetric
	<std::complex<double>, 3,
	 reorder_diagonalize_symmetric, reorder_diagonalize_symmetric, true>,
    Test_diagonalize_symmetric
	<std::complex<double>, 4,
	 reorder_diagonalize_symmetric, reorder_diagonalize_symmetric, true>,
    Test_diagonalize_symmetric
	<std::complex<double>, 6,
	 reorder_diagonalize_symmetric, reorder_diagonalize_symmetric, true>,

    // use Eigen::SelfAdjointEigenSolver
    Test_diagonalize_symmetric
	<double, 6, diagonalize_symmetric, diagonalize_symmetric>,

    Test_diagonalize_symmetric
	<double, 6,
	 reorder_diagonalize_symmetric, reorder_diagonalize_symmetric, true>
> diagonalize_symmetric_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE
(test_diagonalize_symmetric, T, diagonalize_symmetric_tests)
{
    typedef typename T::S S;
    const std::size_t N = T::N;

    Eigen::Matrix<S, N, N> m = Eigen::Matrix<S, N, N>::Random();
    m = ((m + m.transpose())/2).eval();
    Eigen::Array<double, N, 1> s;
    Eigen::Matrix<std::complex<double>, N, N> u;

    T().fxn(m, s, u);
    Eigen::Matrix<std::complex<double>, N, N> diag = u.adjoint() * m * u.conjugate();

    BOOST_CHECK((s >= 0).all());
    for (std::size_t i = 0; i < N; i++)
	for (std::size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-12);

    if (T::check_ascending_order)
	for (std::size_t i = 0; i < N-1; i++)
	    BOOST_CHECK(s[i] <= s[i+1]);

    T().svs(m, s);
    BOOST_CHECK((s >= 0).all());
    for (std::size_t i = 0; i < N; i++)
	for (std::size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 1e-12);
}
#endif // TEST_LINALG2_PART2

template<class S_, int N_,
	 void fxn_(const Eigen::Matrix<S_, N_, N_>&,
		   Eigen::Array<double, N_, 1>&,
		   Eigen::Matrix<S_, N_, N_> *)>
struct Test_diagonalize_hermitian {
    typedef S_ S;
    enum { N = N_ };
    void fxn(const Eigen::Matrix<S_, N_, N_>& m,
	     Eigen::Array<double, N_, 1>& w,
	     Eigen::Matrix<S_, N_, N_> *z)
    { fxn_(m, w, z); }
};

#ifdef TEST_LINALG2_PART3
typedef boost::mpl::list<
    // use Eigen::SelfAdjointEigenSolver
    Test_diagonalize_hermitian<std::complex<double>, 6, hermitian_eigen>,
    Test_diagonalize_hermitian<double	      , 6, hermitian_eigen>
> diagonalize_hermitian_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE
(test_diagonalize_hermitian, T, diagonalize_hermitian_tests)
{
    typedef typename T::S S;
    const std::size_t N = T::N;

    Eigen::Matrix<S, N, N> m = Eigen::Matrix<S, N, N>::Random();
    m = ((m + m.adjoint())/2).eval();
    Eigen::Array<double, N, 1> w;
    Eigen::Matrix<S, N, N> z;

    T().fxn(m, w, &z);		// following LAPACK convention
    Eigen::Matrix<S, N, N> diag = z.adjoint() * m * z;

    for (std::size_t i = 0; i < N; i++)
	for (std::size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? w(i) : 0)), 1e-12);

    for (std::size_t i = 0; i < N-1; i++)
	BOOST_CHECK(w[i] <= w[i+1]);
}
#endif // TEST_LINALG2_PART3

template<class R_, class S_, int M_, int N_ = M_>
struct Test_fs {
    typedef R_ R;
    typedef S_ S;
    enum { M = M_ };
    enum { N = N_ };
};

#ifdef TEST_LINALG2_PART4
typedef boost::mpl::list<
    Test_fs<double, std::complex<double>, 1>,
    Test_fs<double, std::complex<double>, 6>,
    Test_fs<double, std::complex<double>, 4, 6>,
    Test_fs<double, std::complex<double>, 6, 4>,
    Test_fs<double, double	   , 1>,
    Test_fs<double, double	   , 6>,
    Test_fs<double, double	   , 4, 6>,
    Test_fs<double, double	   , 6, 4>,

    Test_fs<long double, std::complex<long double>, 1>,
    Test_fs<long double, std::complex<long double>, 6>,
    Test_fs<long double, std::complex<long double>, 4, 6>,
    Test_fs<long double, std::complex<long double>, 6, 4>,
    Test_fs<long double, long double	     , 1>,
    Test_fs<long double, long double	     , 6>,
    Test_fs<long double, long double	     , 4, 6>,
    Test_fs<long double, long double	     , 6, 4>
> fs_svd_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_fs_svd, T, fs_svd_tests)
{
    typedef typename T::R R;
    typedef typename T::S S;
    const std::size_t M = T::M;
    const std::size_t N = T::N;
    const R eps = numeric_limits<R>::epsilon();

    Eigen::Matrix<S, M, N> m = Eigen::Matrix<S, M, N>::Random();
    Eigen::Array<R, MIN_(M, N), 1> s;
    Eigen::Matrix<S, M, M> u;
    Eigen::Matrix<S, N, N> v;

    gm2calc::fs_svd(m, s, u, v);		// following SARAH convention
    Eigen::Matrix<S, M, N> sigma = u.conjugate() * m * v.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (Eigen::Index i = 0; i < sigma.rows(); i++)
	for (Eigen::Index j = 0; j < sigma.cols(); j++)
	    BOOST_CHECK_SMALL(abs(sigma(i,j) - (i==j ? s(i) : 0)), 50*eps);

    for (Eigen::Index i = 0; i < s.size()-1; i++)
	BOOST_CHECK(s[i] <= s[i+1]);

    gm2calc::fs_svd(m, s);
    BOOST_CHECK((s >= 0).all());
    for (Eigen::Index i = 0; i < sigma.rows(); i++)
	for (Eigen::Index j = 0; j < sigma.cols(); j++)
	    BOOST_CHECK_SMALL(abs(sigma(i,j) - (i==j ? s(i) : 0)), 50*eps);
}
#endif // TEST_LINALG2_PART4

#ifdef TEST_LINALG2_PART5
typedef boost::mpl::list<
    Test_fs<double, double, 1>,
    Test_fs<double, double, 6>,
    Test_fs<double, double, 4, 6>,
    Test_fs<double, double, 6, 4>,

    Test_fs<long double, long double, 1>,
    Test_fs<long double, long double, 6>,
    Test_fs<long double, long double, 4, 6>,
    Test_fs<long double, long double, 6, 4>
> casting_fs_svd_tests;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_casting_fs_svd, T, casting_fs_svd_tests)
{
    typedef typename T::R R;
    const std::size_t M = T::M;
    const std::size_t N = T::N;
    const R eps = numeric_limits<R>::epsilon();

    Eigen::Matrix<R, M, N> m = Eigen::Matrix<R, M, N>::Random();
    Eigen::Array<R, MIN_(M, N), 1> s;
    Eigen::Matrix<std::complex<R>, M, M> u;
    Eigen::Matrix<std::complex<R>, N, N> v;

    gm2calc::fs_svd(m, s, u, v);		// following SARAH convention
    Eigen::Matrix<std::complex<R>, M, N> sigma = u.conjugate() * m * v.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (Eigen::Index i = 0; i < sigma.rows(); i++)
	for (Eigen::Index j = 0; j < sigma.cols(); j++)
	    BOOST_CHECK_SMALL(abs(sigma(i,j) - (i==j ? s(i) : 0)), 50*eps);

    for (Eigen::Index i = 0; i < s.size()-1; i++)
	BOOST_CHECK(s[i] <= s[i+1]);

    gm2calc::fs_svd(m, s);
    BOOST_CHECK((s >= 0).all());
    for (Eigen::Index i = 0; i < sigma.rows(); i++)
	for (Eigen::Index j = 0; j < sigma.cols(); j++)
	    BOOST_CHECK_SMALL(abs(sigma(i,j) - (i==j ? s(i) : 0)), 50*eps);
}
#endif // TEST_LINALG2_PART5

#ifdef TEST_LINALG2_PART6
typedef boost::mpl::list<
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
> fs_diagonalize_symmetric_tests;

template<typename R, typename S, std::size_t N>
void check_fs_diagonalize_symmetric(Eigen::Matrix<S, N, N> m)
{
    const R eps = numeric_limits<R>::epsilon();

    Eigen::Array<R, N, 1> s;
    Eigen::Matrix<std::complex<R>, N, N> u;

    gm2calc::fs_diagonalize_symmetric(m, s, u);
    Eigen::Matrix<std::complex<R>, N, N> diag = u.conjugate() * m * u.adjoint();

    BOOST_CHECK((s >= 0).all());
    for (std::size_t i = 0; i < N; i++)
	for (std::size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 5000*eps);

    for (std::size_t i = 0; i < N-1; i++)
	BOOST_CHECK(s[i] <= s[i+1]);

    gm2calc::fs_diagonalize_symmetric(m, s);
    BOOST_CHECK((s >= 0).all());
    for (std::size_t i = 0; i < N; i++)
	for (std::size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? s(i) : 0)), 5000*eps);
}

template<typename R> R random_angle()
{
    return 2 * M_PI * rand() / RAND_MAX;
}

template<typename S, int N>
struct RandomUnitary;

template<typename R>
struct RandomUnitary<std::complex<R>, 1> {
    Eigen::Matrix<std::complex<R>, 1, 1> operator()() const {
	return Eigen::Matrix<std::complex<R>, 1, 1>(polar(R(1), random_angle<R>()));
    }
};

template<typename R>
struct RandomUnitary<R, 1> {
    Eigen::Matrix<R, 1, 1> operator()() const {
	return Eigen::Matrix<R, 1, 1>(2 * (rand() % 2) - 1);
    }
};

template<typename R>
struct RandomUnitary<std::complex<R>, 2> {
    Eigen::Matrix<std::complex<R>, 2, 2> operator()() const {
	R b = random_angle<R>();
	R c = random_angle<R>();
	R d = random_angle<R>();
	std::complex<R> p = polar(R(1), c);
	std::complex<R> q = polar(R(1), d);
	Eigen::Matrix<std::complex<R>, 2, 2> u;
	u <<      p *cos(b),      q *sin(b),
	    -conj(q)*sin(b), conj(p)*cos(b);
	u *= RandomUnitary<std::complex<R>, 1>()()(0,0);
	return u;
    }
};

template<typename R>
struct RandomUnitary<R, 2> {
    Eigen::Matrix<R, 2, 2> operator()() const {
	R b = random_angle<R>();
	Eigen::Matrix<R, 2, 2> o;
	o << cos(b), sin(b),
	    -sin(b), cos(b);
	o *= RandomUnitary<R, 1>()()(0,0);
	return o;
    }
};

template<typename S, int N>
struct RandomUnitary {
    Eigen::Matrix<S, N, N> operator()() const {
	Eigen::Matrix<S, N, N> u = Eigen::Matrix<S, N, N>::Identity();
	for (int i = 0; i < N-1; i++)
	    for (int j = i+1; j < N; j++) {
		Eigen::Matrix<S, N, N> s = Eigen::Matrix<S, N, N>::Identity();
		Eigen::Matrix<S, 2, 2> r = RandomUnitary<S, 2>()();
		s(i,i) = r(0,0);
		s(i,j) = r(0,1);
		s(j,i) = r(1,0);
		s(j,j) = r(1,1);
		u *= s;
	    }
	return u;
    }
};

BOOST_AUTO_TEST_CASE_TEMPLATE
(test_fs_diagonalize_symmetric, T, fs_diagonalize_symmetric_tests)
{
    typedef typename T::R R;
    typedef typename T::S S;
    const std::size_t N = T::N;

    Eigen::Matrix<S, N, N> m = Eigen::Matrix<S, N, N>::Random();
    m = ((m + m.transpose())/2).eval();
    check_fs_diagonalize_symmetric<R, S, N>(m);

    Eigen::Matrix<S, N, N> u = RandomUnitary<S, N>()();
    Eigen::Array<R, N, 1> s = Eigen::Array<R, N, 1>::Constant(2);
    m = u.transpose() * s.matrix().asDiagonal() * u;
    check_fs_diagonalize_symmetric<R, S, N>(m);
}
#endif // TEST_LINALG2_PART6

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

BOOST_AUTO_TEST_CASE_TEMPLATE
(test_fs_diagonalize_hermitian, P, fs_diagonalize_hermitian_tests)
{
    typedef typename P::first T;
    typedef typename T::R R;
    typedef typename T::S S;
    const std::size_t N = T::N;
    const R eps = numeric_limits<R>::epsilon();

    Eigen::Matrix<S, N, N> m = Eigen::Matrix<S, N, N>::Random();
    m = ((m + m.adjoint())/2).eval();
    Eigen::Array<R, N, 1> w;
    Eigen::Matrix<S, N, N> z;

    gm2calc::fs_diagonalize_hermitian(m, w, z); // following SARAH convention
    Eigen::Matrix<S, N, N> diag = z * m * z.adjoint();

    for (std::size_t i = 0; i < N; i++)
	for (std::size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? w(i) : 0)), 50000*eps);

    for (std::size_t i = 0; i < N-1; i++)
	BOOST_CHECK(abs(w[i]) <= abs(w[i+1]));

    gm2calc::fs_diagonalize_hermitian(m, w);
    for (std::size_t i = 0; i < N; i++)
	for (std::size_t j = 0; j < N; j++)
	    BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j ? w(i) : 0)), 50000*eps);
}
#endif // TEST_LINALG2_PART7

template<int N>
double angle(const Eigen::Matrix<double, 1, N>& u, const Eigen::Matrix<double, 1, N>& v)
{
    // return std::acos(u.dot(v) / (u.norm() * v.norm()));
    Eigen::Matrix<double, 1, N> diff = u - v;
    Eigen::Matrix<double, 1, N> avg  = (u + v) / 2;
    Eigen::Matrix<double, 1, N> n_avg = avg / avg.norm();
    Eigen::Matrix<double, 1, N> diff_perp = diff - diff.dot(n_avg) * n_avg;
    return (diff_perp / avg.norm()).norm();
}

TEST_CASE("test_fs_svd_errbd_easy")
{
    Eigen::Matrix<double, 4, 3> m;
    // example from http://www.netlib.org/lapack/lug/node96.html
    m << 4, 3, 5,
	 2, 5, 8,
	 3, 6, 10,
	 4, 5, 11;
    Eigen::Array<double, 3, 1> s;
    Eigen::Matrix<double, 4, 4> u;
    Eigen::Matrix<double, 3, 3> v;
    double s_errbd;
    Eigen::Array<double, 3, 1> u_errbd;
    Eigen::Array<double, 3, 1> v_errbd;
    gm2calc::fs_svd(m, s, u, v, s_errbd, u_errbd, v_errbd);

    Eigen::Array<double, 3, 1> s_true;
    s_true <<			// from Mathematica
	1.1426562493907868,
	2.3702095896520476,
	21.049381064460057;
    Eigen::Matrix<double, 4, 4> u_true;
    u_true <<			// from Mathematica
	-0.3970292257397899, -0.36237037279930295, -0.3141834916887612, 0.7825242746241264,
	-0.8553177307381059, 0.4111129872735242, 0.2882149672140517, -0.12786642973767423,
	0.3211028895870782, 0.4553747163651497, 0.5708880572549827, 0.603003837531586,
	0.08770580193070293, 0.7016464154456233, -0.7016464154456233, 0.08770580193070293;
    Eigen::Matrix<double, 3, 3> v_true;
    v_true <<			// from Mathematica
	-0.10966642434622631, -0.8536417831240927, 0.5091846241549653,
	-0.9475388908723534, 0.2445224258983454, 0.20586119963989832,
	0.30023878106517676, 0.4598961723449197, 0.8356746885044375;

    Eigen::Array<double, 3, 1> s_error = s - s_true;
    // BOOST_WARN_GE (s_error.abs().minCoeff(), s_errbd / 10);
    CHECK_LE(s_error.abs().maxCoeff(), s_errbd * 10);

    for (std::size_t i = 0; i < 3; i++) {
	double u_error_1 = angle(u.row(i).eval(),   u_true.row(i) .eval());
	double v_error_1 = angle(v.row(i).eval(),   v_true.row(i) .eval());
	double u_error_2 = angle(u.row(i).eval(), (-u_true.row(i)).eval());
	double v_error_2 = angle(v.row(i).eval(), (-v_true.row(i)).eval());
	double u_error, v_error;
	if (u_error_1 + v_error_1 < u_error_2 + v_error_2) {
	    u_error = u_error_1;
	    v_error = v_error_1;
	}
	else {
	    u_error = u_error_2;
	    v_error = v_error_2;
	}
	std::cout << i << ": u_error=" << u_error << " u_errbd=" << u_errbd[i]
	     << " v_error=" << v_error << " v_errbd=" << v_errbd[i] << '\n';
	// BOOST_WARN_GE (u_error, u_errbd[i] / 10);
	CHECK_LE(u_error, u_errbd[i] * 10);
	// BOOST_WARN_GE (v_error, v_errbd[i] / 10);
	CHECK_LE(v_error, v_errbd[i] * 10);
    }
}

TEST_CASE("test_fs_svd_errbd_hard")
{
    Eigen::Matrix<double, 4, 3> m;
    m << 998, -995, -998,
	 999, -996,  996,
	-998,  995,  997,
	-999,  996, -997;
    Eigen::Array<double, 3, 1> s;
    Eigen::Matrix<double, 4, 4> u;
    Eigen::Matrix<double, 3, 3> v;
    double s_errbd;
    Eigen::Array<double, 3, 1> u_errbd;
    Eigen::Array<double, 3, 1> v_errbd;
    gm2calc::fs_svd(m, s, u, v, s_errbd, u_errbd, v_errbd);

    Eigen::Array<double, 3, 1> s_true;
    s_true <<			// from Mathematica
	1.0670512753940286e-6,
	1994.0005015055856,
	2819.945389541342;
    Eigen::Matrix<double, 4, 4> u_true;
    u_true <<			// from Mathematica
	-0.4997490015980897, -0.5002505150161606, -0.5002508726439843, -0.4997493592259134,
	-0.500501377470506, 0.4994983710239794, 0.49999987308309773, -0.49999987541138774,
	-0.49974918603411456, -0.5002506882142286, 0.49974918603354523, 0.5002506882136593,
	-0.5, 0.5, -0.5, 0.5;
    Eigen::Matrix<double, 3, 3> v_true;
    v_true <<			// from Mathematica
	0.7060421306428913, 0.7081698311535926, 1.0670511412087742e-6,
	-7.521989593869741e-7, -7.567975556452352e-7, 0.9999999999994309,
	-0.708169831153997, 0.7060421306432921, 1.605386026477888e-9;

    Eigen::Array<double, 3, 1> s_error = s - s_true;
    // BOOST_WARN_GE (s_error.abs().minCoeff(), s_errbd / 10);
    CHECK_LE(s_error.abs().maxCoeff(), s_errbd * 10);

    for (std::size_t i = 0; i < 3; i++) {
	double u_error_1 = angle(u.row(i).eval(),   u_true.row(i) .eval());
	double v_error_1 = angle(v.row(i).eval(),   v_true.row(i) .eval());
	double u_error_2 = angle(u.row(i).eval(), (-u_true.row(i)).eval());
	double v_error_2 = angle(v.row(i).eval(), (-v_true.row(i)).eval());
	double u_error, v_error;
	if (u_error_1 + v_error_1 < u_error_2 + v_error_2) {
	    u_error = u_error_1;
	    v_error = v_error_1;
	}
	else {
	    u_error = u_error_2;
	    v_error = v_error_2;
	}
	std::cout << i << ": u_error=" << u_error << " u_errbd=" << u_errbd[i]
	     << " v_error=" << v_error << " v_errbd=" << v_errbd[i] << '\n';
	// this m seems to be a very bad matrix for error estimation
	// for singular vectors
	// BOOST_WARN_GE (u_error, u_errbd[i] / 10);
	CHECK_LE(u_error, u_errbd[i] * 1e5);
	// BOOST_WARN_GE (v_error, v_errbd[i] / 10);
	CHECK_LE(v_error, v_errbd[i] * 1e5);
    }
}

TEST_CASE("test_fs_diagonalize_hermitian_errbd")
{
    Eigen::Matrix<double, 3, 3> m;
    // example from http://www.netlib.org/lapack/lug/node89.html
    m << 1, 2, 3,
	 2, 4, 5,
	 3, 5, 6;
    Eigen::Array<double, 3, 1> w;
    Eigen::Matrix<double, 3, 3> z;
    double w_errbd;
    Eigen::Array<double, 3, 1> z_errbd;
    gm2calc::fs_diagonalize_hermitian(m, w, z, w_errbd, z_errbd);

    Eigen::Array<double, 3, 1> w_true;
    w_true <<			// from Mathematica
	0.1709151888271795,
	-0.5157294715892571,
	11.34481428276208;
    Eigen::Matrix<double, 3, 3> z_true;
    z_true <<			// from Mathematica
	0.5910090485061035, -0.7369762290995782, 0.3279852776056818,
	-0.7369762290995782, -0.3279852776056818, 0.5910090485061035,
	0.3279852776056818, 0.5910090485061035, 0.7369762290995782;

    Eigen::Array<double, 3, 1> w_error = w - w_true;
    // BOOST_WARN_GE (w_error.abs().minCoeff(), w_errbd / 10);
    CHECK_LE(w_error.abs().maxCoeff(), w_errbd * 10);

    for (std::size_t i = 0; i < 3; i++) {
	double z_error_1 = angle(z.row(i).eval(),   z_true.row(i) .eval());
	double z_error_2 = angle(z.row(i).eval(), (-z_true.row(i)).eval());
	double z_error = std::min(z_error_1, z_error_2);
	std::cout << i << ": z_error=" << z_error << " z_errbd=" << z_errbd[i]
	     << '\n';
	// BOOST_WARN_GE (z_error, z_errbd[i] / 10);
	CHECK_LE(z_error, z_errbd[i] * 10);
    }
}

TEST_CASE("test_diagonalize_symmetric_errbd")
{
   const std::complex<double> m14(-4.4299349683288838e-06,7.0893778912230837e-06);
   const double m24 = 4.1620179585196247e-05;
   const double m34 = 3.7871449517912927e-05;
   const double m44 = 90.229999999964136;
   Eigen::Matrix<std::complex<double>, 4, 4> m;
   m <<
      0. , 0. , 0. , m14,
      0. , 0. , 0. , m24,
      0. , 0. , 0. , m34,
      m14, m24, m34, m44;

   Eigen::Array<double, 4, 1> s;
   Eigen::Matrix<std::complex<double>, 4, 4> u;
   double s_errbd = 0;
   Eigen::Array<double, 4, 1> u_errbd;
   gm2calc::fs_diagonalize_symmetric_errbd(m, s, &u, &s_errbd, &u_errbd);
   INFO("u =\n" << u);
   INFO("uu^+ =\n" << (u*u.adjoint()).eval());
   CHECK(u.isUnitary());
}
