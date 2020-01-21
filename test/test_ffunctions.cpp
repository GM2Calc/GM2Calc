#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "ffunctions.hpp"
#include <limits>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)

TEST_CASE("limits -> 0")
{
   const double PI = 3.141592653589793;
   const double eps = std::numeric_limits<double>::epsilon();

   using namespace gm2calc;

   CHECK_CLOSE(F1C(0.0), 4.0, eps);
   CHECK_CLOSE(F1N(0.0), 2.0, eps);
   CHECK_CLOSE(F2N(0.0), 3.0, eps);
   CHECK_CLOSE(F3N(0.0), 8.0/105.0, eps);
   CHECK_CLOSE(F4N(0.0),  -3.0*(-9.0 + PI*PI)/4.0, eps);
}

TEST_CASE("limits -> 1")
{
   const double eps = std::numeric_limits<double>::epsilon();

   using namespace gm2calc;

   CHECK_CLOSE(F1C(1.0), 1.0, eps);
   CHECK_CLOSE(F1N(1.0), 1.0, eps);
   CHECK_CLOSE(F2N(1.0), 1.0, eps);
   CHECK_CLOSE(F3N(1.0), 1.0, eps);
   CHECK_CLOSE(F4N(1.0), 1.0, eps);

   CHECK_CLOSE(Fa(1.0, 1.0), 0.25, eps);
   CHECK_CLOSE(Fb(1.0, 1.0), 1.0/12.0, eps);
}
