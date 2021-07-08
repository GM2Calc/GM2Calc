// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2calc/GeneralTHDM.hpp"

#include <iostream>

namespace gm2calc {

/// Table 1, arxiv:1607.06292
double GeneralTHDM::get_zeta_u() const
{
   return 1.0/get_tan_beta();
}

/// Table 1, arxiv:1607.06292
double GeneralTHDM::get_zeta_d() const
{
   switch (yukawa_scheme) {
      case Yukawa_scheme::type_1:
         return 1.0/get_tan_beta();
      case Yukawa_scheme::type_2:
         return -get_tan_beta();
      case Yukawa_scheme::type_X:
         return 1.0/get_tan_beta();
      case Yukawa_scheme::type_Y:
         return -get_tan_beta();
   }
}

/// Table 1, arxiv:1607.06292
double GeneralTHDM::get_zeta_l() const
{
   switch (yukawa_scheme) {
      case Yukawa_scheme::type_1:
         return 1.0/get_tan_beta();
      case Yukawa_scheme::type_2:
         return -get_tan_beta();
      case Yukawa_scheme::type_X:
         return -get_tan_beta();
      case Yukawa_scheme::type_Y:
         return 1.0/get_tan_beta();
   }
}

std::ostream& operator<<(std::ostream& ostr, const GeneralTHDM& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace gm2calc
