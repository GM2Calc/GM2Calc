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

#include "GeneralTHDM_soft_parameters.hpp"

#include <iostream>

namespace gm2calc {

void GeneralTHDM_soft_parameters::print(std::ostream& ostr) const
{
   GeneralTHDM_susy_parameters::print(ostr);
   ostr << "dimensionfull parameters at Q = " << get_scale() << ":\n";
   ostr << "M122 = " << M122 << '\n';
   ostr << "M112 = " << M112 << '\n';
   ostr << "M222 = " << M222 << '\n';
   ostr << "v1 = " << v1 << '\n';
   ostr << "v2 = " << v2 << '\n';
}

std::ostream& operator<<(std::ostream& ostr, const GeneralTHDM_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace gm2calc
