#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_ffunctions.hpp"
#include "read_data.hpp"
#include <cmath>
#include <limits>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)

TEST_CASE("limits -> 0")
{
   const double PI = 3.141592653589793;
   const double eps = std::numeric_limits<double>::epsilon();

   CHECK_CLOSE(gm2calc::F1C(0.0), 4.0, eps);
   CHECK_CLOSE(gm2calc::F1N(0.0), 2.0, eps);
   CHECK_CLOSE(gm2calc::F2N(0.0), 3.0, eps);
   CHECK_CLOSE(gm2calc::F3N(0.0), 8.0/105.0, eps);
   CHECK_CLOSE(gm2calc::F4N(0.0), -3.0*(-9.0 + PI*PI)/4.0, eps);

   CHECK(gm2calc::Iabc(0.0, 0.0, 0.0) == 0.0);
   CHECK(gm2calc::Iabc(0.0, 0.0, 1.0) == 0.0);
   CHECK_CLOSE(gm2calc::Iabc(0.0, 1.0, 1.0), 1.0, eps);
   CHECK_CLOSE(gm2calc::Iabc(0.0, 1.0, 2.0), std::log(4.0)/3, eps);
}

TEST_CASE("limits -> 1")
{
   const double eps = std::numeric_limits<double>::epsilon();

   CHECK_CLOSE(gm2calc::F1C(1.0), 1.0, eps);
   CHECK_CLOSE(gm2calc::F2C(1.0), 1.0, eps);
   CHECK_CLOSE(gm2calc::F3C(1.0), 1.0, eps);
   CHECK_CLOSE(gm2calc::F4C(1.0), 1.0, eps);

   CHECK_CLOSE(gm2calc::F1N(1.0), 1.0, eps);
   CHECK_CLOSE(gm2calc::F2N(1.0), 1.0, eps);
   CHECK_CLOSE(gm2calc::F3N(1.0), 1.0, eps);
   CHECK_CLOSE(gm2calc::F4N(1.0), 1.0, eps);

   CHECK_CLOSE(gm2calc::Fa(1.0, 1.0), 0.25, eps);
   CHECK_CLOSE(gm2calc::Fb(1.0, 1.0), 1.0/12.0, eps);

   CHECK_CLOSE(gm2calc::G3(1.0), 1.0/3.0, eps);
   CHECK_CLOSE(gm2calc::G4(1.0), 1.0/6.0, eps);

   CHECK_CLOSE(gm2calc::Iabc(1.0, 1.0, 1.0), 0.5, eps);
}

TEST_CASE("symmetries")
{
   const double eps = std::numeric_limits<double>::epsilon();

   using gm2calc::Iabc;

   CHECK_CLOSE(gm2calc::Fa(2.0, 1.0), gm2calc::Fa(1.0, 2.0), eps);
   CHECK_CLOSE(gm2calc::Fa(2.0, 3.0), gm2calc::Fa(3.0, 2.0), eps);
   CHECK_CLOSE(gm2calc::Fb(2.0, 1.0), gm2calc::Fb(1.0, 2.0), eps);
   CHECK_CLOSE(gm2calc::Fb(2.0, 3.0), gm2calc::Fb(3.0, 2.0), eps);

   CHECK_CLOSE(Iabc(1.0, 1.0, 2.0), Iabc(1.0, 2.0, 1.0), eps);
   CHECK_CLOSE(Iabc(1.0, 1.0, 2.0), Iabc(2.0, 1.0, 1.0), eps);
   CHECK_CLOSE(Iabc(1.0, 2.0, 2.0), Iabc(2.0, 2.0, 1.0), eps);
   CHECK_CLOSE(Iabc(1.0, 2.0, 2.0), Iabc(2.0, 1.0, 2.0), eps);
   CHECK_CLOSE(Iabc(1.0, 2.0, 3.0), Iabc(2.0, 3.0, 1.0), eps);
   CHECK_CLOSE(Iabc(1.0, 2.0, 3.0), Iabc(3.0, 2.0, 1.0), eps);
   CHECK_CLOSE(Iabc(1.0, 2.0, 3.0), Iabc(2.0, 1.0, 3.0), eps);
   CHECK_CLOSE(Iabc(1.0, 2.0, 3.0), Iabc(1.0, 3.0, 2.0), eps);
   CHECK_CLOSE(Iabc(1.0, 2.0, 3.0), Iabc(3.0, 2.0, 1.0), eps);
}

template <typename T>
void test_1(const char* func_name, T func, double eps)
{
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + func_name + ".txt");
   const auto data = gm2calc::test::read_from_file<double>(filename);

   for (const auto& v: data) {
      if (v.size() < 2) {
         continue;
      }

      const auto x          = v.at(0);
      const auto f_expected = v.at(1);
      const auto f_gm2calc  = func(x);

      INFO("x = " << x << ", " << func_name << "(expected) = " << f_expected
                  << ", " << func_name << "(gm2calc) = " << f_gm2calc);

      CHECK_CLOSE(f_expected, f_gm2calc, eps);
   }
}

TEST_CASE("test_F")
{
   test_1("F1C", [] (double x) { return gm2calc::F1C(x); }, 1e-9);
   test_1("F2C", [] (double x) { return gm2calc::F2C(x); }, 1e-9);
   test_1("F3C", [] (double x) { return gm2calc::F3C(x); }, 1e-9);
   test_1("F4C", [] (double x) { return gm2calc::F4C(x); }, 1e-9);

   test_1("F1N", [] (double x) { return gm2calc::F1N(x); }, 1e-9);
   test_1("F2N", [] (double x) { return gm2calc::F2N(x); }, 1e-9);
   test_1("F3N", [] (double x) { return gm2calc::F3N(x); }, 1e-9);
   test_1("F4N", [] (double x) { return gm2calc::F4N(x); }, 1e-9);

   test_1("F1" , [] (double x) { return gm2calc::F1(x);  }, 1e-14);
   test_1("F1t", [] (double x) { return gm2calc::F1t(x); }, 1e-14);
   test_1("F2" , [] (double x) { return gm2calc::F2(x);  }, 1e-14);
   test_1("F3" , [] (double x) { return gm2calc::F3(x);  }, 1e-14);

   test_1("fPS", [] (double x) { return gm2calc::f_PS(x);}, 1e-14);
   test_1("fS" , [] (double x) { return gm2calc::f_S(x); }, 1e-14);
   test_1("fsferm", [] (double x) { return gm2calc::f_sferm(x); }, 1e-14);
}

TEST_CASE("test_G")
{
   test_1("G3", [] (double x) { return gm2calc::G3(x); }, 1e-9);
   test_1("G4", [] (double x) { return gm2calc::G4(x); }, 1e-9);
}

template <typename T>
void test_2(const char* func_name, T func, double eps)
{
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + func_name + ".txt");
   const auto data = gm2calc::test::read_from_file<double>(filename);

   for (const auto& v: data) {
      if (v.size() < 3) {
         continue;
      }

      const auto x          = v.at(0);
      const auto y          = v.at(1);
      const auto f_expected = v.at(2);
      const auto f_gm2calc  = func(x, y);

      INFO("x = " << x << ", y = " << y << ", " << func_name
                  << "(expected) = " << f_expected << ", " << func_name
                  << "(gm2calc) = " << f_gm2calc);

      CHECK_CLOSE(f_expected, f_gm2calc, eps);
   }
}

TEST_CASE("test_Fab")
{
   test_2("Fa", [] (double x, double y) { return gm2calc::Fa(x,y); }, 1e-7);
   test_2("Fb", [] (double x, double y) { return gm2calc::Fb(x,y); }, 1e-7);
}

template <typename T>
void test_3(const char* func_name, T func, double eps)
{
   const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR + func_name + ".txt");
   const auto data = gm2calc::test::read_from_file<double>(filename);

   for (const auto& v: data) {
      if (v.size() < 4) {
         continue;
      }

      const auto x          = v.at(0);
      const auto y          = v.at(1);
      const auto z          = v.at(2);
      const auto f_expected = v.at(3);
      const auto f_gm2calc  = func(x, y, z);

      INFO("x = " << x << ", y = " << y << ", z = " << z << ", " << func_name
                  << "(expected) = " << f_expected << ", " << func_name
                  << "(gm2calc) = " << f_gm2calc);

      CHECK_CLOSE(f_expected, f_gm2calc, eps);
   }
}

TEST_CASE("test_Iabc")
{
   test_3("Iabc", [] (double x, double y, double z) { return gm2calc::Iabc(x,y,z); }, 1e-7);
   test_3("Phi" , [] (double x, double y, double z) { return gm2calc::Phi(x,y,z);  }, 1e-14);
}
