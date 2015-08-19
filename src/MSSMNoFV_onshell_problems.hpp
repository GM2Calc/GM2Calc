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

#ifndef MSSMNoFV_onshell_PROBLEMS_H
#define MSSMNoFV_onshell_PROBLEMS_H

#include <iosfwd>
#include <string>

namespace gm2calc {

class MSSMNoFV_onshell_problems {
public:
   MSSMNoFV_onshell_problems();
   void clear();
   void clear_problems();
   void clear_warnings();
   void flag_no_convergence_Mu_MassB_MassWB(bool flag = true);
   void flag_tachyon(const std::string&);
   bool have_problem() const;
   bool have_warning() const;
   std::string get_warning() const;
   void print(std::ostream&) const;
   void print_problems(std::ostream&) const;
   void print_warnings(std::ostream&) const;

private:
   bool have_no_convergence_Mu_MassB_MassWB;
   bool have_tachyon;
   std::string tachyonic_particle;
};

std::ostream& operator<<(std::ostream&, const MSSMNoFV_onshell_problems&);

} // namespace gm2calc

#endif
