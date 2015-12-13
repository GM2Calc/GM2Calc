#ifndef TEST_H
#define TEST_H

#include <iostream>
#include <string>
#include <cmath>
#include <limits>

namespace gm2calc_test {

static unsigned long g_error = 0;
static bool g_verbose = false;

#define VERBOSE_MSG(msg)                          \
   do {                                           \
      if (g_verbose) {                            \
         std::cout << msg << std::endl;           \
      }                                           \
   } while (0)

#define CHECK_EQUAL(a, b) \
   check_equal((a), (b), #a, #b)

#define CHECK_CLOSE_ABS(a, b, eps) \
   check_close_abs((a), (b), (eps), #a, #b)

#define CHECK_CLOSE_REL(a, b, eps) \
   check_close_rel((a), (b), (eps), #a, #b)

template <typename T>
bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon())
{
   return std::abs(a) < prec;
}

template <typename T>
bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
{
   return is_zero(a - b, prec);
}

template <typename T>
bool is_equal(long a, long b, T prec = std::numeric_limits<T>::epsilon())
{
   return is_zero(a - b, prec);
}

template <typename T>
bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon())
{
   if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
      return true;

   if (std::abs(a) < std::numeric_limits<T>::epsilon())
      return is_equal(a, b, prec);

   return std::abs((a - b)/a) < prec;
}

void check_equal(double a, double b, const std::string& astr, const std::string& bstr)
{
   if (a == b) {
      VERBOSE_MSG("test passed: " << astr << " == " << bstr
                  << ": " << a << " == " << b);
   } else {
      std::cout << "test failed: " << astr << " == " << bstr
                << ": " << a << " != " << b << std::endl;
      g_error++;
   }
}

void check_close_abs(double a, double b, double eps, const std::string& astr, const std::string& bstr)
{
   if (is_equal(a, b, eps)) {
      VERBOSE_MSG("test passed: " << astr << " ~ " << bstr
                  << ": " << a << " ~ " << b
                  << " (required maximum absolute difference: " << eps << ")");
   } else {
      std::cout << "test failed: " << astr << " ~ " << bstr
                << ": " << a << " !~ " << b
                << " (required maximum absolute difference: " << eps << ")"
                << std::endl;
      g_error++;
   }
}

void check_close_rel(double a, double b, double eps, const std::string& astr, const std::string& bstr)
{
   if (is_equal_rel(a, b, eps)) {
      VERBOSE_MSG("test passed: " << astr << " ~ " << bstr
                  << ": " << a << " ~ " << b
                  << " (required maximum relative difference: " << eps);
   } else {
      std::cout << "test failed: " << astr << " ~ " << bstr
                << ": " << a << " !~ " << b
                << " (required maximum relative difference: " << eps << ")"
                << std::endl;
      g_error++;
   }
}

} // namespace gm2calc_test

#endif
