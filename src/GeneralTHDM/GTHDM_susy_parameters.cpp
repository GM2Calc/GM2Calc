// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "GTHDM_susy_parameters.hpp"

#include <iostream>

namespace gm2calc {

void GTHDM_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "dimensionless parameters at Q = " << get_scale() << ":\n";
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "Lambda6 = " << Lambda6 << '\n';
   ostr << "Lambda5 = " << Lambda5 << '\n';
   ostr << "Lambda7 = " << Lambda7 << '\n';
   ostr << "Lambda1 = " << Lambda1 << '\n';
   ostr << "Lambda4 = " << Lambda4 << '\n';
   ostr << "Lambda3 = " << Lambda3 << '\n';
   ostr << "Lambda2 = " << Lambda2 << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "Xu = " << Xu << '\n';
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "Xd = " << Xd << '\n';
   ostr << "Xe = " << Xe << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const GTHDM_susy_parameters& susy_pars)
{
   susy_pars.print(ostr);
   return ostr;
}

} // namespace gm2calc
