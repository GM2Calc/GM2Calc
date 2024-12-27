#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_numerics.hpp"


TEST_CASE("sign")
{
   CHECK_EQ(gm2calc::sign(2), 1.0);
   CHECK_EQ(gm2calc::sign(-2), -1.0);
   CHECK_EQ(gm2calc::sign(0), 1.0);
}


TEST_CASE("signed_sqr")
{
   CHECK_EQ(gm2calc::signed_sqr(2), 4.0);
   CHECK_EQ(gm2calc::signed_sqr(-2), -4.0);
}


TEST_CASE("signed_abs_sqrt")
{
   CHECK_EQ(gm2calc::signed_abs_sqrt(4), 2.0);
   CHECK_EQ(gm2calc::signed_abs_sqrt(-4), -2.0);
}
