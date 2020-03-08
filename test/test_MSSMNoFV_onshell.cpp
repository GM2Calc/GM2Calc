#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_linalg.hpp"
#include <Eigen/Core>


namespace {

template <class Derived>
unsigned find_right_like_smuon(const Eigen::MatrixBase<Derived>& ZM)
{
   // Incorrect:
   // return (ZM(0,0) > ZM(0,1)) ? 1 : 0;

   return (std::abs(ZM(0,0)) > std::abs(ZM(0,1))) ? 1 : 0;
}

} // anonymous namespace


TEST_CASE("right_like_smuon")
{
   Eigen::Matrix<double,2,2> M;
   Eigen::Matrix<double,2,2> Z;
   Eigen::Array<double,2,1> m;

   M(0,0) = 1.0;
   M(0,1) = 0.1;
   M(1,0) = M(0,1);
   M(1,1) = 2.0;

   gm2calc::fs_diagonalize_hermitian(M, m, Z);

   const auto idx = find_right_like_smuon(Z);
   const auto mR = m(idx); // right-like smuoon

   INFO("Z = " << Z);
   INFO("m = " << m.transpose());
   INFO("Mass of right-like smuon [" << idx << "] = " << mR);

   CHECK(std::abs(mR - M(1,1)) <= std::abs(mR - M(0,0)));
}
