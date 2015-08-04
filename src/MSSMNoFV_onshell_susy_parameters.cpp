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

#include "MSSMNoFV_onshell_susy_parameters.hpp"

#include <iostream>

namespace flexiblesusy {
namespace gm2calc {

MSSMNoFV_onshell_susy_parameters::MSSMNoFV_onshell_susy_parameters()
   : scale(0)
   , Yd(Eigen::Matrix<double,3,3>::Zero()), Ye(Eigen::Matrix<double,3,3>::Zero(
   )), Yu(Eigen::Matrix<double,3,3>::Zero()), Mu(0), g1(0), g2(0), g3(0), vd(0)
   , vu(0)
{
}

MSSMNoFV_onshell_susy_parameters::MSSMNoFV_onshell_susy_parameters(
   double scale_
   , const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_
   , const Eigen::Matrix<double,3,3>& Yu_, double Mu_, double g1_, double g2_,
   double g3_, double vd_, double vu_

)
   : scale(scale_)
   , Yd(Yd_), Ye(Ye_), Yu(Yu_), Mu(Mu_), g1(g1_), g2(g2_), g3(g3_), vd(vd_), vu
   (vu_)
{
}

void MSSMNoFV_onshell_susy_parameters::clear()
{
   scale = 0.;
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   Yu = Eigen::Matrix<double,3,3>::Zero();
   Mu = 0.;
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   vd = 0.;
   vu = 0.;
}

void MSSMNoFV_onshell_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters:\n";
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "Mu = " << Mu << '\n';
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "vd = " << vd << '\n';
   ostr << "vu = " << vu << '\n';
}

std::ostream& operator<<(std::ostream& ostr, const MSSMNoFV_onshell_susy_parameters& susy_pars)
{
   susy_pars.print(std::cout);
   return ostr;
}

} // namespace gm2calc
} // namespace flexiblesusy
