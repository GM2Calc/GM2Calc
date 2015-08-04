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

#include "MSSMNoFV_onshell_problems.hpp"
#include "numerics2.hpp"

#include <iostream>

namespace flexiblesusy {
namespace gm2calc {

MSSMNoFV_onshell_problems::MSSMNoFV_onshell_problems()
   : have_tachyon(false)
   , tachyonic_particle("")
{
}

void MSSMNoFV_onshell_problems::clear()
{
   have_tachyon = false;
   tachyonic_particle.clear();
}

void MSSMNoFV_onshell_problems::flag_tachyon(const std::string& particle_name)
{
   have_tachyon = true;
   tachyonic_particle = particle_name;
}

bool MSSMNoFV_onshell_problems::have_problem() const
{
   return have_tachyon;
}

void MSSMNoFV_onshell_problems::print(std::ostream& ostr) const
{
   if (have_tachyon)
      ostr << "Problem: " << tachyonic_particle << " tachyon";
}

std::ostream& operator<<(std::ostream& ostr, const MSSMNoFV_onshell_problems& problems)
{
   problems.print(ostr);
   return ostr;
}

} // namespace gm2calc
} // namespace flexiblesusy
