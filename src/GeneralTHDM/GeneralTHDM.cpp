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
#include "gm2calc/gm2_error.hpp"
#include "gm2_numerics.hpp"

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

void GeneralTHDM::set_basis(const GeneralTHDM::General_basis& basis)
{
   set_Lambda1(basis.lambda1);
   set_Lambda2(basis.lambda2);
   set_Lambda3(basis.lambda3);
   set_Lambda4(basis.lambda4);
   set_Lambda5(basis.lambda5);
   set_Lambda6(basis.lambda6);
   set_Lambda7(basis.lambda7);
   set_tan_beta_and_v(basis.tan_beta, get_v_from_mW(sm.get_mw()));
   set_M122(basis.M122);

   solve_ewsb();
   calculate_MSbar_masses();
}

void GeneralTHDM::set_basis(const GeneralTHDM::Physical_basis& basis)
{
   if (basis.mh > basis.mH) {
      throw EInvalidInput("mh must be less than or equal to mH.");
   }

   if (basis.tan_beta <= 0) {
      throw EInvalidInput("tan(beta) must be greater than zero.");
   }

   if (std::abs(basis.sin_beta_minus_alpha) > 1) {
      throw EInvalidInput("|sin(beta - alpha_h)| must be less than or equal to 1.");
   }

   if (std::abs(basis.mh) < 0) {
      throw EInvalidInput("mh must be greater than or equal to zero.");
   }

   if (std::abs(basis.mH) < 0) {
      throw EInvalidInput("mH must be greater than or equal to zero.");
   }

   if (std::abs(basis.mA) < 0) {
      throw EInvalidInput("mA must be greater than or equal to zero.");
   }

   if (std::abs(basis.mHp) < 0) {
      throw EInvalidInput("mHp must be greater than or equal to zero.");
   }

   const double sba = basis.sin_beta_minus_alpha;
   const double tb = basis.tan_beta;
   const double ctb = 1./tb;
   const double sb = tb/std::sqrt(1 + sqr(tb));
   const double cb = 1./std::sqrt(1 + sqr(tb));
   const double sb2 = sqr(sb);
   const double cb2 = sqr(cb);
   const double alpha = -std::asin(sba) + std::atan(tb);
   const double sa = std::sin(alpha);
   const double ca  = std::cos(alpha);
   const double sa2 = sqr(sa);
   const double ca2 = sqr(ca);
   const double mh = basis.mh;
   const double mH = basis.mH;
   const double mA = basis.mA;
   const double mHp = basis.mHp;
   const double lambda6 = basis.lambda6;
   const double lambda7 = basis.lambda7;
   const double m12_2 = basis.M122;
   const double v = get_v_from_mW(sm.get_mw());
   const double v2 = sqr(v);

   set_Lambda1((sqr(mH)*ca2 + sqr(mh)*sa2 - m12_2*tb)/v2/cb2 - 1.5*lambda6*tb + 0.5*lambda7*tb*tb*tb);
   set_Lambda2((sqr(mH)*sa2 + sqr(mh)*ca2 - m12_2*ctb)/v2/sb2 + 0.5*lambda6*ctb*ctb*ctb - 1.5*lambda7*ctb);
   set_Lambda3(((sqr(mH) - sqr(mh))*ca*sa + 2.*mHp*mHp*sb*cb - m12_2)/v2/sb/cb - 0.5*lambda6*ctb - 0.5*lambda7*tb);
   set_Lambda4(((sqr(mA) - 2.*mHp*mHp)*cb*sb + m12_2)/v2/sb/cb - 0.5*lambda6*ctb - 0.5*lambda7*tb);
   set_Lambda5((m12_2 - sqr(mA)*sb*cb)/v2/sb/cb - 0.5*lambda6*ctb - 0.5*lambda7*tb);
   set_Lambda6(basis.lambda6);
   set_Lambda7(basis.lambda7);
   set_tan_beta_and_v(basis.tan_beta, v);
   set_M122(basis.M122);

   solve_ewsb();
   calculate_MSbar_masses();
}

void GeneralTHDM::set_sm(const SM& sm_)
{
   sm = sm_;
   set_alpha_em_and_cw(sm.get_alpha_em_mz(), sm.get_cw());
}

void GeneralTHDM::set_tan_beta(double tb)
{
   set_tan_beta_and_v(tb, sm.get_v());
}

std::ostream& operator<<(std::ostream& ostr, const GeneralTHDM& model)
{
   model.print(ostr);

   ostr << "Yukawa scheme: ";

   switch (model.get_yukawa_scheme()) {
      case GeneralTHDM::Yukawa_scheme::type_1:
         ostr << "Type I\n";
         break;
      case GeneralTHDM::Yukawa_scheme::type_2:
         ostr << "Type II\n";
         break;
      case GeneralTHDM::Yukawa_scheme::type_X:
         ostr << "Type X\n";
         break;
      case GeneralTHDM::Yukawa_scheme::type_Y:
         ostr << "Type Y\n";
         break;
   }

   ostr << "zeta_u = " << model.get_zeta_u() << '\n';
   ostr << "zeta_d = " << model.get_zeta_d() << '\n';
   ostr << "zeta_l = " << model.get_zeta_l() << '\n';

   return ostr;
}

} // namespace gm2calc
