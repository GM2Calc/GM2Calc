#ifndef TEST_H
#define TEST_H

#include <iostream>
#include <string>

namespace gm2calc_test {

static int g_error = 0;
static bool g_verbose = false;

#define VERBOSE_MSG(msg)                          \
   do {                                           \
      if (g_verbose) {                            \
         std::cout << msg << std::endl;           \
      }                                           \
   } while (0)

#define CHECK_EQUAL(a, b) \
   check_equal((a), (b), #a, #b)

void check_equal(double a, double b, const std::string& astr, const std::string& bstr)
{
   if (a == b) {
      VERBOSE_MSG("test passed: " << astr << " == " << bstr
                  << ": " << a << " == " << b);
   } else {
      std::cout << "test failed: " << astr << " == " << bstr
                << ": " << a << " != " << b << std::endl;
      g_error = 1;
   }
}

} // namespace gm2calc_test

#endif
