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

#include "gm2calc/GeneralTHDM_parameters.hpp"

#include <iostream>

namespace gm2calc {

void GeneralTHDM_parameters::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "Parameters:\n"
           "----------------------------------------\n";
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "lambda1 = " << lambda1 << '\n';
   ostr << "lambda2 = " << lambda2 << '\n';
   ostr << "lambda3 = " << lambda3 << '\n';
   ostr << "lambda4 = " << lambda4 << '\n';
   ostr << "lambda5 = " << lambda5 << '\n';
   ostr << "lambda6 = " << lambda6 << '\n';
   ostr << "lambda7 = " << lambda7 << '\n';
   ostr << "Yu =\n" << Yu << '\n';
   ostr << "Xu =\n" << Xu << '\n';
   ostr << "Yd =\n" << Yd << '\n';
   ostr << "Yl =\n" << Yl << '\n';
   ostr << "Xd =\n" << Xd << '\n';
   ostr << "Xl =\n" << Xl << '\n';
   ostr << "m122 = " << m122 << '\n';
   ostr << "m112 = " << m112 << '\n';
   ostr << "m222 = " << m222 << '\n';
   ostr << "v1 = " << v1 << '\n';
   ostr << "v2 = " << v2 << '\n';
}

std::ostream& operator<<(std::ostream& ostr, const GeneralTHDM_parameters& pars)
{
   pars.print(ostr);
   return ostr;
}

} // namespace gm2calc
