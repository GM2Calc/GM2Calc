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

#include "gm2calc/GeneralTHDM_problems.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>

namespace gm2calc {

void GeneralTHDM_problems::clear_problems()
{
   tachyons.clear();
}

void GeneralTHDM_problems::clear_warnings()
{
}

void GeneralTHDM_problems::clear()
{
   clear_problems();
   clear_warnings();
}

void GeneralTHDM_problems::flag_tachyon(const std::string& particle_name)
{
   tachyons.push_back(particle_name);
   std::sort(tachyons.begin(), tachyons.end());
   tachyons.erase(std::unique(tachyons.begin(), tachyons.end()), tachyons.end());
}

bool GeneralTHDM_problems::have_tachyon() const
{
   return !tachyons.empty();
}

bool GeneralTHDM_problems::have_problem() const
{
   return have_tachyon();
}

bool GeneralTHDM_problems::have_warning() const
{
   return false;
}

std::string GeneralTHDM_problems::get_warnings() const
{
   std::ostringstream ostr;
   print_warnings(ostr);
   return ostr.str();
}

std::string GeneralTHDM_problems::get_problems() const
{
   std::ostringstream ostr;
   print_problems(ostr);
   return ostr.str();
}

void GeneralTHDM_problems::print_problems(std::ostream& ostr) const
{
   if (have_problem()) {
      ostr << "Problem: ";
   }

   if (have_tachyon()) {
      for (auto it = tachyons.cbegin(), end = tachyons.cend(); it != end; ++it) {
         if (it != tachyons.begin()) {
            ostr << ", ";
         }
         ostr << *it << " tachyon";
      }
   }
}

void GeneralTHDM_problems::print_warnings(std::ostream& ostr) const
{
   if (have_warning()) {
      ostr << "Warning:";
   }
}

void GeneralTHDM_problems::print(std::ostream& ostr) const
{
   print_warnings(ostr);
   print_problems(ostr);
}

std::ostream& operator<<(std::ostream& ostr, const GeneralTHDM_problems& problems)
{
   problems.print(ostr);
   return ostr;
}

} // namespace gm2calc
