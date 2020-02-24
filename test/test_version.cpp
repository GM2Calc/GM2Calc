#include "gm2calc/gm2_version.h"

#include <iostream>

int main()
{
   std::cout << "GM2Calc " << GM2CALC_VERSION << '\n'
             << "major: " << GM2CALC_VERSION_MAJOR << '\n'
             << "minor: " << GM2CALC_VERSION_MINOR << '\n'
             << "patch: " << GM2CALC_VERSION_PATCH << std::endl;

   return 0;
}
