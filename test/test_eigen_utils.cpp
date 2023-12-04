#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_eigen_utils.hpp"


TEST_CASE("normalize_to_interval")
{
   Eigen::Matrix<double,2,2> m;
   m << 2.0, 3.0, -2.0, -3.0;

   gm2calc::normalize_to_interval(m, -1.0, 1.0);

   CHECK_EQ(m(0,0), 1.0);
   CHECK_EQ(m(0,1), 1.0);
   CHECK_EQ(m(1,0), -1.0);
   CHECK_EQ(m(1,1), -1.0);
}


TEST_CASE("remove_if_equal")
{
   Eigen::Array<double,4,1> src;
   src << 1.0, 2.0, 3.0, 4.0;

   Eigen::Array<double,2,1> cmp;
   cmp << 1.1, 2.9;

   const Eigen::Array<double,2,1> res = gm2calc::remove_if_equal(src, cmp);

   CHECK_EQ(res(0), 2.0);
   CHECK_EQ(res(1), 4.0);
}
