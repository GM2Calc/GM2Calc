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

#include "gm2calc/MSSMNoFV_onshell_soft_parameters.hpp"

#include <iostream>
#include <utility>

namespace gm2calc {

MSSMNoFV_onshell_soft_parameters::MSSMNoFV_onshell_soft_parameters()
   : TYd(Eigen::Matrix<double,3,3>::Zero()), TYe(Eigen::Matrix<double,3,3>
   ::Zero()), TYu(Eigen::Matrix<double,3,3>::Zero()), BMu(0), mq2(Eigen::Matrix
   <double,3,3>::Zero()), ml2(Eigen::Matrix<double,3,3>::Zero()), mHd2(0), mHu2
   (0), md2(Eigen::Matrix<double,3,3>::Zero()), mu2(Eigen::Matrix<double,3,3>
   ::Zero()), me2(Eigen::Matrix<double,3,3>::Zero()), MassB(0), MassWB(0),
   MassG(0)
{
}

MSSMNoFV_onshell_soft_parameters::MSSMNoFV_onshell_soft_parameters(
   const MSSMNoFV_onshell_susy_parameters& susy_model
   , Eigen::Matrix<double,3,3>  TYd_, Eigen::Matrix<double,3,3> 
   TYe_, Eigen::Matrix<double,3,3>  TYu_, double BMu_, Eigen::Matrix<double,3,3>  mq2_, Eigen::Matrix<double,3,3>  ml2_,
   double mHd2_, double mHu2_, Eigen::Matrix<double,3,3>  md2_, Eigen::Matrix<double,3,3>  mu2_, Eigen::Matrix<double,3,3>  me2_,
   double MassB_, double MassWB_, double MassG_

)
   : MSSMNoFV_onshell_susy_parameters(susy_model)
   , TYd(std::move(TYd_)), TYe(std::move(TYe_)), TYu(std::move(TYu_)), BMu(BMu_), mq2(std::move(mq2_)), ml2(std::move(ml2_)), mHd2(
   mHd2_), mHu2(mHu2_), md2(std::move(md2_)), mu2(std::move(mu2_)), me2(std::move(me2_)), MassB(MassB_), MassWB(
   MassWB_), MassG(MassG_)
{
}

void MSSMNoFV_onshell_soft_parameters::clear()
{
   MSSMNoFV_onshell_susy_parameters::clear();

   TYd = Eigen::Matrix<double,3,3>::Zero();
   TYe = Eigen::Matrix<double,3,3>::Zero();
   TYu = Eigen::Matrix<double,3,3>::Zero();
   BMu = 0.;
   mq2 = Eigen::Matrix<double,3,3>::Zero();
   ml2 = Eigen::Matrix<double,3,3>::Zero();
   mHd2 = 0.;
   mHu2 = 0.;
   md2 = Eigen::Matrix<double,3,3>::Zero();
   mu2 = Eigen::Matrix<double,3,3>::Zero();
   me2 = Eigen::Matrix<double,3,3>::Zero();
   MassB = 0.;
   MassWB = 0.;
   MassG = 0.;
}

void MSSMNoFV_onshell_soft_parameters::print(std::ostream& ostr) const
{
   MSSMNoFV_onshell_susy_parameters::print(ostr);
   ostr << "soft parameters:\n";
   ostr << "TYd = " << TYd << '\n';
   ostr << "TYe = " << TYe << '\n';
   ostr << "TYu = " << TYu << '\n';
   ostr << "BMu = " << BMu << '\n';
   ostr << "mq2 = " << mq2 << '\n';
   ostr << "ml2 = " << ml2 << '\n';
   ostr << "mHd2 = " << mHd2 << '\n';
   ostr << "mHu2 = " << mHu2 << '\n';
   ostr << "md2 = " << md2 << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "me2 = " << me2 << '\n';
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "MassG = " << MassG << '\n';
}

std::ostream& operator<<(std::ostream& ostr, const MSSMNoFV_onshell_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace gm2calc
