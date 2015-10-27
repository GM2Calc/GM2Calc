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
#include <sstream>
#include <algorithm>

namespace gm2calc {

MSSMNoFV_onshell_problems::MSSMNoFV_onshell_problems()
   : have_no_convergence_Mu_MassB_MassWB(false)
   , have_no_convergence_me2(false)
   , tachyons()
   , convergence_problem_Mu_MassB_MassWB()
   , convergence_problem_me2()
{
}

void MSSMNoFV_onshell_problems::clear_problems()
{
   tachyons.clear();
}

void MSSMNoFV_onshell_problems::clear_warnings()
{
   have_no_convergence_Mu_MassB_MassWB = false;
   have_no_convergence_me2 = false;
   convergence_problem_Mu_MassB_MassWB.clear();
   convergence_problem_me2.clear();
}

void MSSMNoFV_onshell_problems::clear()
{
   clear_problems();
   clear_warnings();
}

void MSSMNoFV_onshell_problems::flag_tachyon(const std::string& particle_name)
{
   tachyons.push_back(particle_name);
   std::sort(tachyons.begin(), tachyons.end());
   tachyons.erase(std::unique(tachyons.begin(), tachyons.end()), tachyons.end());
}

void MSSMNoFV_onshell_problems::flag_no_convergence_Mu_MassB_MassWB(
   double precision, unsigned iterations)
{
   have_no_convergence_Mu_MassB_MassWB = true;
   convergence_problem_Mu_MassB_MassWB.precision = precision;
   convergence_problem_Mu_MassB_MassWB.iterations = iterations;
}

void MSSMNoFV_onshell_problems::unflag_no_convergence_Mu_MassB_MassWB()
{
   have_no_convergence_Mu_MassB_MassWB = false;
   convergence_problem_Mu_MassB_MassWB.clear();
}

void MSSMNoFV_onshell_problems::flag_no_convergence_me2(
   double precision, unsigned iterations)
{
   have_no_convergence_me2 = true;
   convergence_problem_me2.precision = precision;
   convergence_problem_me2.iterations = iterations;
}

void MSSMNoFV_onshell_problems::unflag_no_convergence_me2()
{
   have_no_convergence_me2 = false;
   convergence_problem_me2.clear();
}

bool MSSMNoFV_onshell_problems::have_tachyon() const
{
   return !tachyons.empty();
}

bool MSSMNoFV_onshell_problems::have_problem() const
{
   return have_tachyon();
}

bool MSSMNoFV_onshell_problems::have_warning() const
{
   return have_no_convergence_Mu_MassB_MassWB || have_no_convergence_me2;
}

std::string MSSMNoFV_onshell_problems::get_warnings() const
{
   std::ostringstream ostr;
   print_warnings(ostr);
   return ostr.str();
}

std::string MSSMNoFV_onshell_problems::get_problems() const
{
   std::ostringstream ostr;
   print_problems(ostr);
   return ostr.str();
}

void MSSMNoFV_onshell_problems::print_problems(std::ostream& ostr) const
{
   if (have_problem())
      ostr << "Problem: ";

   if (have_tachyon()) {
      for (std::vector<std::string>::const_iterator it = tachyons.begin(),
              end = tachyons.end(); it != end; ++it) {
         if (it != tachyons.begin())
            ostr << ", ";
         ostr << *it << " tachyon";
      }
   }
}

void MSSMNoFV_onshell_problems::print_warnings(std::ostream& ostr) const
{
   if (have_warning())
      ostr << "Warning:";

   if (have_no_convergence_Mu_MassB_MassWB)
      ostr << " DR-bar to on-shell conversion for Mu, M1, M2 failed"
              " (reached absolute accuracy: "
           << convergence_problem_Mu_MassB_MassWB.precision << " GeV),";

   if (have_no_convergence_me2)
      ostr << " DR-bar to on-shell conversion for me2 failed"
              " (reached absolute accuracy: "
           << convergence_problem_me2.precision << " GeV)";
}

void MSSMNoFV_onshell_problems::print(std::ostream& ostr) const
{
   print_warnings(ostr);
   print_problems(ostr);
}

std::ostream& operator<<(std::ostream& ostr, const MSSMNoFV_onshell_problems& problems)
{
   problems.print(ostr);
   return ostr;
}

} // namespace gm2calc
