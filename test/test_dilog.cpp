#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_dilog.hpp"


#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


TEST_CASE("real dilog")
{
   const struct Data {
      double x;
      double y;
   } data[] = {
      {-5.0, -2.749279126060808290025587515376268644497},
      {-4.0, -2.369939796998365831985537425350323048751},
      {-3.0, -1.939375420766708953077271719177891441223},
      {-2.0, -1.43674636688368094636290202389358335425 },
      {-1.0, -0.822467033424113218236207583323012594609},
      {-0.5, -0.448414206923646202443064405915774320834},
      { 0.0,  0.0},
      { 0.2,  0.211003775439704772611185096074072517246},
      { 0.5,  0.582240526465012505902656320159680108744},
      { 0.7,  0.889377624286038738601006274807361793537},
      { 1.0,  1.644934066848226436472415166646025189219},
      { 1.2,  2.129169430383959659444305569819388370407},
      { 2.2,  2.458586601999741851737986678758721928694}
   };

   for (const auto& d: data) {
      CHECK_CLOSE(gm2calc::dilog(d.x), d.y, 1e-15);
   }
}


TEST_CASE("complex dilog")
{
   const double eps = 1e-15;

   CHECK_EQ(gm2calc::dilog(0.0), 0.0);
}


TEST_CASE("real clausen_2")
{
   const double eps = 1e-15;

   CHECK_EQ(gm2calc::clausen_2(0.0), 0.0);
}
