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

#include "gm2calc/THDM.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2_numerics.hpp"

namespace gm2calc {

namespace {

const double sqrt2 = 1.4142135623730950; // Sqrt[2]

/// Eq.(6) arxiv:0908.1554 and arxiv:1001.0293, solved for xi_f
double calc_xi_bar(double zeta, double tan_beta) noexcept
{
   return (tan_beta + zeta)/(1 - tan_beta*zeta);
}

} // anonymous namespace

THDM::THDM(const thdm::Gauge_basis& basis, const SM& sm_)
   : sm(sm_)
   , yukawa_scheme(basis.yukawa_scheme)
   , zeta_u(basis.zeta_u)
   , zeta_d(basis.zeta_d)
   , zeta_l(basis.zeta_l)
{
   init_gauge_couplings();
   set_basis(basis);
}

THDM::THDM(const thdm::Mass_basis& basis, const SM& sm_)
   : sm(sm_)
   , yukawa_scheme(basis.yukawa_scheme)
   , zeta_u(basis.zeta_u)
   , zeta_d(basis.zeta_d)
   , zeta_l(basis.zeta_l)
{
   init_gauge_couplings();
   set_basis(basis);
}

void THDM::init_gauge_couplings()
{
   set_alpha_em_and_cw(sm.get_alpha_em_mz(), sm.get_cw());
}

/**
 * Initialize Yukawa couplings
 */
void THDM::init_yukawas()
{
   const Eigen::Matrix<double,3,3> mu = sm.get_mu().asDiagonal();
   const Eigen::Matrix<double,3,3> md = sm.get_md().asDiagonal();
   const Eigen::Matrix<double,3,3> ml = sm.get_ml().asDiagonal();

   switch (yukawa_scheme) {
   case thdm::Yukawa_scheme::type_1:
      Xu.setZero();
      Xd.setZero();
      Xl.setZero();
      set_Yu(sqrt2*mu/v1);
      set_Yd(sqrt2*md/v1);
      set_Yl(sqrt2*ml/v1);
      break;
   case thdm::Yukawa_scheme::type_2:
      Yu.setZero();
      Xd.setZero();
      Xl.setZero();
      set_Xu(sqrt2*mu/v2);
      set_Yd(sqrt2*md/v1);
      set_Yl(sqrt2*ml/v1);
      break;
   case thdm::Yukawa_scheme::type_X:
      Yu.setZero();
      Yd.setZero();
      Xl.setZero();
      set_Xu(sqrt2*mu/v2);
      set_Xd(sqrt2*md/v2);
      set_Yl(sqrt2*ml/v1);
      break;
   case thdm::Yukawa_scheme::type_Y:
      Xu.setZero();
      Xd.setZero();
      Yl.setZero();
      set_Yu(sqrt2*mu/v1);
      set_Yd(sqrt2*md/v1);
      set_Xl(sqrt2*ml/v2);
      break;
   case thdm::Yukawa_scheme::aligned:
      set_Yu(sqrt2*mu/(v1 + v2*calc_xi_bar(get_zeta_u(), get_tan_beta())));
      set_Yd(sqrt2*md/(v1 + v2*calc_xi_bar(get_zeta_d(), get_tan_beta())));
      set_Yl(sqrt2*ml/(v1 + v2*calc_xi_bar(get_zeta_l(), get_tan_beta())));
      set_Xu(calc_xi_bar(get_zeta_u(), get_tan_beta())*Yu);
      set_Xd(calc_xi_bar(get_zeta_d(), get_tan_beta())*Yd);
      set_Xl(calc_xi_bar(get_zeta_l(), get_tan_beta())*Yl);
      break;
   case thdm::Yukawa_scheme::general:
      set_Yu(sqrt2*mu/v1 - v2/v1*Xu);
      set_Yd(sqrt2*md/v1 - v2/v1*Xd);
      set_Yl(sqrt2*ml/v1 - v2/v1*Xl);
      break;
   }
}

/// Table 1, arxiv:1607.06292
double THDM::get_zeta_u() const
{
   switch (yukawa_scheme) {
   case thdm::Yukawa_scheme::aligned:
      return zeta_u;
   default:
      return 1.0/get_tan_beta();
   }
}

/// Table 1, arxiv:1607.06292
double THDM::get_zeta_d() const
{
   switch (yukawa_scheme) {
   case thdm::Yukawa_scheme::type_1:
      return 1.0/get_tan_beta();
   case thdm::Yukawa_scheme::type_2:
      return -get_tan_beta();
   case thdm::Yukawa_scheme::type_X:
      return 1.0/get_tan_beta();
   case thdm::Yukawa_scheme::type_Y:
      return -get_tan_beta();
   case thdm::Yukawa_scheme::aligned:
      return zeta_d;
   case thdm::Yukawa_scheme::general:
      return -get_tan_beta(); // should never arrive here
   }
}

/// Table 1, arxiv:1607.06292
double THDM::get_zeta_l() const
{
   switch (yukawa_scheme) {
   case thdm::Yukawa_scheme::type_1:
      return 1.0/get_tan_beta();
   case thdm::Yukawa_scheme::type_2:
      return -get_tan_beta();
   case thdm::Yukawa_scheme::type_X:
      return -get_tan_beta();
   case thdm::Yukawa_scheme::type_Y:
      return 1.0/get_tan_beta();
   case thdm::Yukawa_scheme::aligned:
      return zeta_l;
   case thdm::Yukawa_scheme::general:
      return -get_tan_beta(); // should never arrive here
   }
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_xi_u() const
{
   const double cb = get_cos_beta();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> mu = sm.get_mu().asDiagonal();

   if (yukawa_scheme != thdm::Yukawa_scheme::general) {
      return sqrt2*mu*get_zeta_u()/v;
   }

   return get_Xu().real()/cb - sqrt2*mu*get_tan_beta()/v;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_xi_d() const
{
   const double cb = get_cos_beta();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> md = sm.get_md().asDiagonal();

   if (yukawa_scheme != thdm::Yukawa_scheme::general) {
      return sqrt2*md*get_zeta_d()/v;
   }

   return get_Xd()/cb - sqrt2*md*get_tan_beta()/v;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_xi_l() const
{
   const double cb = get_cos_beta();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> ml = sm.get_ml().asDiagonal();

   if (yukawa_scheme != thdm::Yukawa_scheme::general) {
      return sqrt2*ml*get_zeta_l()/v;
   }

   return get_Xl()/cb - sqrt2*ml*get_tan_beta()/v;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_yuh() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> mu = sm.get_mu().asDiagonal();

   return sba*mu/v + cba*get_xi_u()/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_yuH() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> mu = sm.get_mu().asDiagonal();

   return cba*mu/v - sba*get_xi_u()/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_yuA() const
{
   return get_xi_u()/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_yuHp() const
{
   return -(sm.get_ckm().adjoint()*get_xi_u()).transpose();
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ydh() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> md = sm.get_md().asDiagonal();

   return sba*md/v + cba*get_xi_d()/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ydH() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> md = sm.get_md().asDiagonal();

   return cba*md/v - sba*get_xi_d()/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ydA() const
{
   return -get_xi_d()/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ydHp() const
{
   return sm.get_ckm()*get_xi_d();
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ylh() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> ml = sm.get_ml().asDiagonal();

   return sba*ml/v + cba*get_xi_l()/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ylH() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> ml = sm.get_ml().asDiagonal();

   return cba*ml/v - sba*get_xi_l()/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ylA() const
{
   return -get_xi_l()/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ylHp() const
{
   return get_xi_l();
}

void THDM::set_basis(const thdm::Gauge_basis& basis)
{
   set_lambda1(basis.lambda1);
   set_lambda2(basis.lambda2);
   set_lambda3(basis.lambda3);
   set_lambda4(basis.lambda4);
   set_lambda5(basis.lambda5);
   set_lambda6(basis.lambda6);
   set_lambda7(basis.lambda7);
   set_tan_beta_and_v(basis.tan_beta, sm.get_v());
   set_m122(basis.m122);
   set_Xu(basis.Xu);
   set_Xd(basis.Xd);
   set_Xl(basis.Xl);

   init_yukawas();
   solve_ewsb();
   calculate_MSbar_masses();

   if (get_problems().have_problem()) {
      throw EPhysicalProblem(get_problems().get_problems());
   }
}

void THDM::set_basis(const thdm::Mass_basis& basis)
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
   const double m12_2 = basis.m122;
   const double v = sm.get_v();
   const double v2 = sqr(v);

   set_lambda1((sqr(mH)*ca2 + sqr(mh)*sa2 - m12_2*tb)/v2/cb2 - 1.5*lambda6*tb + 0.5*lambda7*tb*tb*tb);
   set_lambda2((sqr(mH)*sa2 + sqr(mh)*ca2 - m12_2*ctb)/v2/sb2 + 0.5*lambda6*ctb*ctb*ctb - 1.5*lambda7*ctb);
   set_lambda3(((sqr(mH) - sqr(mh))*ca*sa + 2.*mHp*mHp*sb*cb - m12_2)/v2/sb/cb - 0.5*lambda6*ctb - 0.5*lambda7*tb);
   set_lambda4(((sqr(mA) - 2.*mHp*mHp)*cb*sb + m12_2)/v2/sb/cb - 0.5*lambda6*ctb - 0.5*lambda7*tb);
   set_lambda5((m12_2 - sqr(mA)*sb*cb)/v2/sb/cb - 0.5*lambda6*ctb - 0.5*lambda7*tb);
   set_lambda6(basis.lambda6);
   set_lambda7(basis.lambda7);
   set_tan_beta_and_v(basis.tan_beta, v);
   set_m122(basis.m122);
   set_Xu(basis.Xu);
   set_Xd(basis.Xd);
   set_Xl(basis.Xl);

   init_yukawas();
   solve_ewsb();
   calculate_MSbar_masses();

   if (get_problems().have_problem()) {
      throw EPhysicalProblem(get_problems().get_problems());
   }
}

void THDM::set_tan_beta(double tb)
{
   set_tan_beta_and_v(tb, sm.get_v());
}

std::ostream& operator<<(std::ostream& ostr, const THDM& model)
{
   model.print(ostr);

   ostr << "Yukawa scheme: ";

   switch (model.get_yukawa_scheme()) {
      case thdm::Yukawa_scheme::type_1:
         ostr << "Type I\n";
         break;
      case thdm::Yukawa_scheme::type_2:
         ostr << "Type II\n";
         break;
      case thdm::Yukawa_scheme::type_X:
         ostr << "Type X\n";
         break;
      case thdm::Yukawa_scheme::type_Y:
         ostr << "Type Y\n";
         break;
      case thdm::Yukawa_scheme::aligned:
         ostr << "Aligned\n";
         break;
      case thdm::Yukawa_scheme::general:
         ostr << "General\n";
         break;
   }

   ostr << "zeta_u = " << model.get_zeta_u() << '\n';
   ostr << "zeta_d = " << model.get_zeta_d() << '\n';
   ostr << "zeta_l = " << model.get_zeta_l() << '\n';

   return ostr;
}

} // namespace gm2calc
