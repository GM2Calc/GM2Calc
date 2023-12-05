#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_numerics.hpp"


TEST_CASE("sign")
{
   CHECK_EQ(gm2calc::sign(2), 1.0);
   CHECK_EQ(gm2calc::sign(-2), -1.0);
   CHECK_EQ(gm2calc::sign(0), 1.0);
}
