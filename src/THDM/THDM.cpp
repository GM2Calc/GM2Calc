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
#include "gm2_log.hpp"
#include "gm2_mf.hpp"
#include "gm2_numerics.hpp"

#include <cmath>

namespace gm2calc {

namespace {

const double sqrt2 = 1.4142135623730950; // Sqrt[2]

/// Eq.(6) arxiv:0908.1554 and arxiv:1001.0293, solved for xi_f
double calc_xi_bar(double zeta, double tan_beta) noexcept
{
   return (tan_beta + zeta)/(1 - tan_beta*zeta);
}

} // anonymous namespace

THDM::THDM(const thdm::Gauge_basis& basis, const SM& sm_, const thdm::Config& cfg)
   : sm(sm_)
   , yukawa_type(basis.yukawa_type)
   , zeta_u(basis.zeta_u)
   , zeta_d(basis.zeta_d)
   , zeta_l(basis.zeta_l)
   , config(cfg)
{
   init_gauge_couplings();
   set_basis(basis);
}

THDM::THDM(const thdm::Mass_basis& basis, const SM& sm_, const thdm::Config& cfg)
   : sm(sm_)
   , yukawa_type(basis.yukawa_type)
   , zeta_u(basis.zeta_u)
   , zeta_d(basis.zeta_d)
   , zeta_l(basis.zeta_l)
   , config(cfg)
{
   init_gauge_couplings();
   set_basis(basis);
}

/// returns the up-type quark masses
Eigen::Matrix<double,3,1> THDM::get_mu(double scale) const
{
   Eigen::Matrix<double,3,1> mu = sm.get_mu();

   if (config.running_couplings && scale > 0) {
      const double mt_pole = mu(2);
      // replace mt_pole by mt(SM(6), MS-bar, Q = scale)
      mu(2) = calculate_mt_SM6_MSbar(mt_pole, sm.get_alpha_s_mz(), sm.get_mz(), scale);
   }

   return mu;
}

/// returns the down-type quark masses
Eigen::Matrix<double,3,1> THDM::get_md(double scale) const
{
   Eigen::Matrix<double,3,1> md = sm.get_md();

   if (config.running_couplings && scale > 0) {
      const double mb_mb = md(2);
      const double mt_pole = sm.get_mu(2);
      // replace mb_mb by mb(SM(6), MS-bar, Q = scale)
      md(2) = calculate_mb_SM6_MSbar(mb_mb, mt_pole, sm.get_alpha_s_mz(), sm.get_mz(), scale);
   }

   return md;
}

/// returns the charged lepton masses
Eigen::Matrix<double,3,1> THDM::get_ml(double scale) const
{
   Eigen::Matrix<double,3,1> ml = sm.get_ml();

   if (config.running_couplings && scale > 0) {
      const double mtau_pole = ml(2);
      // replace mtau_pole by mtau(SM(6), MS-bar, Q = scale)
      ml(2) = calculate_mtau_SM6_MSbar(mtau_pole, sm.get_alpha_em_mz(), scale);
   }

   return ml;
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
   const Eigen::Matrix<std::complex<double>,3,3> vckm_adj = sm.get_ckm().adjoint();

   switch (yukawa_type) {
   case thdm::Yukawa_type::type_1:
      Xu.setZero();
      Xd.setZero();
      Xl.setZero();
      set_Yu(sqrt2*vckm_adj*mu/v1);
      set_Yd(sqrt2*md/v1);
      set_Yl(sqrt2*ml/v1);
      break;
   case thdm::Yukawa_type::type_2:
      Yu.setZero();
      Xd.setZero();
      Xl.setZero();
      set_Xu(sqrt2*vckm_adj*mu/v2);
      set_Yd(sqrt2*md/v1);
      set_Yl(sqrt2*ml/v1);
      break;
   case thdm::Yukawa_type::type_X:
      Yu.setZero();
      Yd.setZero();
      Xl.setZero();
      set_Xu(sqrt2*vckm_adj*mu/v2);
      set_Xd(sqrt2*md/v2);
      set_Yl(sqrt2*ml/v1);
      break;
   case thdm::Yukawa_type::type_Y:
      Xu.setZero();
      Xd.setZero();
      Yl.setZero();
      set_Yu(sqrt2*vckm_adj*mu/v1);
      set_Yd(sqrt2*md/v1);
      set_Xl(sqrt2*ml/v2);
      break;
   case thdm::Yukawa_type::aligned:
      set_Yu(sqrt2*vckm_adj*mu/(v1 + v2*calc_xi_bar(get_zeta_u(), get_tan_beta())));
      set_Yd(sqrt2*md/(v1 + v2*calc_xi_bar(get_zeta_d(), get_tan_beta())));
      set_Yl(sqrt2*ml/(v1 + v2*calc_xi_bar(get_zeta_l(), get_tan_beta())));
      set_Xu(calc_xi_bar(get_zeta_u(), get_tan_beta())*Yu);
      set_Xd(calc_xi_bar(get_zeta_d(), get_tan_beta())*Yd);
      set_Xl(calc_xi_bar(get_zeta_l(), get_tan_beta())*Yl);
      break;
   case thdm::Yukawa_type::general:
      set_Yu(sqrt2*vckm_adj*mu/v1 - v2/v1*Xu);
      set_Yd(sqrt2*md/v1 - v2/v1*Xd);
      set_Yl(sqrt2*ml/v1 - v2/v1*Xl);
      break;
   }
}

/// Table 1, arxiv:1607.06292
double THDM::get_zeta_u() const
{
   switch (yukawa_type) {
   case thdm::Yukawa_type::type_1:
   case thdm::Yukawa_type::type_2:
   case thdm::Yukawa_type::type_X:
   case thdm::Yukawa_type::type_Y:
      return 1.0/get_tan_beta();
   case thdm::Yukawa_type::aligned:
      return zeta_u;
   case thdm::Yukawa_type::general:
      return 0; // never used in the general THDM
   }
   throw ESetupError("Bug: unhandled case in get_zeta_u.");
}

/// Table 1, arxiv:1607.06292
double THDM::get_zeta_d() const
{
   switch (yukawa_type) {
   case thdm::Yukawa_type::type_1:
      return 1.0/get_tan_beta();
   case thdm::Yukawa_type::type_2:
      return -get_tan_beta();
   case thdm::Yukawa_type::type_X:
      return 1.0/get_tan_beta();
   case thdm::Yukawa_type::type_Y:
      return -get_tan_beta();
   case thdm::Yukawa_type::aligned:
      return zeta_d;
   case thdm::Yukawa_type::general:
      return 0; // never used in the general THDM
   }
   throw ESetupError("Bug: unhandled case in get_zeta_d.");
}

/// Table 1, arxiv:1607.06292
double THDM::get_zeta_l() const
{
   switch (yukawa_type) {
   case thdm::Yukawa_type::type_1:
      return 1.0/get_tan_beta();
   case thdm::Yukawa_type::type_2:
      return -get_tan_beta();
   case thdm::Yukawa_type::type_X:
      return -get_tan_beta();
   case thdm::Yukawa_type::type_Y:
      return 1.0/get_tan_beta();
   case thdm::Yukawa_type::aligned:
      return zeta_l;
   case thdm::Yukawa_type::general:
      return 0; // never used in the general THDM
   }
   throw ESetupError("Bug: unhandled case in get_zeta_l.");
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_xi_u(const Eigen::Matrix<double,3,3>& mu) const
{
   const double cb = get_cos_beta();
   const double v = get_v();

   if (yukawa_type != thdm::Yukawa_type::general) {
      return sqrt2*mu*get_zeta_u()/v;
   }

   return get_Xu().real()/cb - sqrt2*mu*get_tan_beta()/v;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_xi_d(const Eigen::Matrix<double,3,3>& md) const
{
   const double cb = get_cos_beta();
   const double v = get_v();

   if (yukawa_type != thdm::Yukawa_type::general) {
      return sqrt2*md*get_zeta_d()/v;
   }

   return get_Xd()/cb - sqrt2*md*get_tan_beta()/v;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_xi_l(const Eigen::Matrix<double,3,3>& ml) const
{
   const double cb = get_cos_beta();
   const double v = get_v();

   if (yukawa_type != thdm::Yukawa_type::general) {
      return sqrt2*ml*get_zeta_l()/v;
   }

   return get_Xl()/cb - sqrt2*ml*get_tan_beta()/v;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_yuh() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> mu = get_mu(get_Mhh(0)).asDiagonal();

   return sba*mu/v + cba*get_xi_u(mu)/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_yuH() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> mu = get_mu(get_Mhh(1)).asDiagonal();

   return cba*mu/v - sba*get_xi_u(mu)/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_yuA() const
{
   const Eigen::Matrix<double,3,3> mu = get_mu(get_MAh(1)).asDiagonal();
   return get_xi_u(mu)/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_yuHp() const
{
   const Eigen::Matrix<double,3,3> mu = get_mu(get_MHm(1)).asDiagonal();
   return -(sm.get_ckm().adjoint()*get_xi_u(mu)).transpose();
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ydh() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> md = get_md(get_Mhh(0)).asDiagonal();

   return sba*md/v + cba*get_xi_d(md)/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ydH() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> md = get_md(get_Mhh(1)).asDiagonal();

   return cba*md/v - sba*get_xi_d(md)/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ydA() const
{
   const Eigen::Matrix<double,3,3> md = get_md(get_MAh(1)).asDiagonal();
   return -get_xi_d(md)/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ydHp() const
{
   const Eigen::Matrix<double,3,3> md = get_md(get_MHm(1)).asDiagonal();
   return sm.get_ckm()*get_xi_d(md);
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ylh() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> ml = get_ml(get_Mhh(0)).asDiagonal();

   return sba*ml/v + cba*get_xi_l(ml)/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ylH() const
{
   const double cba = get_cos_beta_minus_alpha();
   const double sba = get_sin_beta_minus_alpha();
   const double v = get_v();
   const Eigen::Matrix<double,3,3> ml = get_ml(get_Mhh(1)).asDiagonal();

   return cba*ml/v - sba*get_xi_l(ml)/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ylA() const
{
   const Eigen::Matrix<double,3,3> ml = get_ml(get_MAh(1)).asDiagonal();
   return -get_xi_l(ml)/sqrt2;
}

Eigen::Matrix<std::complex<double>,3,3> THDM::get_ylHp() const
{
   const Eigen::Matrix<double,3,3> ml = get_ml(get_MHm(1)).asDiagonal();
   return get_xi_l(ml);
}

void THDM::set_basis(const thdm::Gauge_basis& basis)
{
   if (basis.tan_beta <= 0) {
      const char* err = "tan(beta) must be greater than zero.";
      if (config.force_output) {
         WARNING(err);
      } else {
         throw EInvalidInput(err);
      }
   }

   set_lambda1(basis.lambda(0));
   set_lambda2(basis.lambda(1));
   set_lambda3(basis.lambda(2));
   set_lambda4(basis.lambda(3));
   set_lambda5(basis.lambda(4));
   set_lambda6(basis.lambda(5));
   set_lambda7(basis.lambda(6));
   set_tan_beta_and_v(basis.tan_beta, sm.get_v());
   set_m122(basis.m122);
   set_Xu(basis.Xu);
   set_Xd(basis.Xd);
   set_Xl(basis.Xl);

   validate();

   solve_ewsb();
   init_yukawas();
   calculate_MSbar_masses();

   if (get_problems().have_problem()) {
      if (config.force_output) {
         WARNING(get_problems().get_problems());
      } else {
         throw EPhysicalProblem(get_problems().get_problems());
      }
   }
}

void THDM::set_basis(const thdm::Mass_basis& basis)
{
   if (basis.mh > basis.mH) {
      const char* err = "mh must be less than or equal to mH.";
      if (config.force_output) {
         WARNING(err);
      } else {
         throw EInvalidInput(err);
      }
   }

   if (basis.tan_beta <= 0) {
      const char* err = "tan(beta) must be greater than zero.";
      if (config.force_output) {
         WARNING(err);
      } else {
         throw EInvalidInput(err);
      }
   }

   if (std::abs(basis.sin_beta_minus_alpha) > 1) {
      const char* err = "|sin(beta - alpha_h)| must be less than or equal to 1.";
      if (config.force_output) {
         WARNING(err);
      } else {
         throw EInvalidInput(err);
      }
   }

   if (basis.mh < 0) {
      const char* err = "mh must be greater than or equal to zero.";
      if (config.force_output) {
         WARNING(err);
      } else {
         throw EInvalidInput(err);
      }
   }

   if (basis.mH < 0) {
      const char* err = "mH must be greater than or equal to zero.";
      if (config.force_output) {
         WARNING(err);
      } else {
         throw EInvalidInput(err);
      }
   }

   if (basis.mA < 0) {
      const char* err = "mA must be greater than or equal to zero.";
      if (config.force_output) {
         WARNING(err);
      } else {
         throw EInvalidInput(err);
      }
   }

   if (basis.mHp < 0) {
      const char* err = "mHp must be greater than or equal to zero.";
      if (config.force_output) {
         WARNING(err);
      } else {
         throw EInvalidInput(err);
      }
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
   const double lambda6 = basis.lambda_6;
   const double lambda7 = basis.lambda_7;
   const double m12_2 = basis.m122;
   const double v = sm.get_v();
   const double v2 = sqr(v);

   set_lambda1((sqr(mH)*ca2 + sqr(mh)*sa2 - m12_2*tb)/v2/cb2 - 1.5*lambda6*tb + 0.5*lambda7*tb*tb*tb);
   set_lambda2((sqr(mH)*sa2 + sqr(mh)*ca2 - m12_2*ctb)/v2/sb2 + 0.5*lambda6*ctb*ctb*ctb - 1.5*lambda7*ctb);
   set_lambda3(((sqr(mH) - sqr(mh))*ca*sa + 2.*mHp*mHp*sb*cb - m12_2)/v2/sb/cb - 0.5*lambda6*ctb - 0.5*lambda7*tb);
   set_lambda4(((sqr(mA) - 2.*mHp*mHp)*cb*sb + m12_2)/v2/sb/cb - 0.5*lambda6*ctb - 0.5*lambda7*tb);
   set_lambda5((m12_2 - sqr(mA)*sb*cb)/v2/sb/cb - 0.5*lambda6*ctb - 0.5*lambda7*tb);
   set_lambda6(basis.lambda_6);
   set_lambda7(basis.lambda_7);
   set_tan_beta_and_v(basis.tan_beta, v);
   set_m122(basis.m122);
   set_Xu(basis.Xu);
   set_Xd(basis.Xd);
   set_Xl(basis.Xl);

   validate();

   solve_ewsb();
   init_yukawas();
   calculate_MSbar_masses();

   if (get_problems().have_problem()) {
      if (config.force_output) {
         WARNING(get_problems().get_problems());
      } else {
         throw EPhysicalProblem(get_problems().get_problems());
      }
   }
}

void THDM::set_tan_beta(double tb)
{
   set_tan_beta_and_v(tb, sm.get_v());
}

void THDM::print(std::ostream& ostr) const
{
   THDM_mass_eigenstates::print(ostr);

   ostr << "Yukawa type: " << yukawa_type_to_string() << '\n';
   ostr << "Running couplings: " << (config.running_couplings ? "yes" : "no") << '\n';

   ostr << "zeta_u = " << get_zeta_u() << '\n';
   ostr << "zeta_d = " << get_zeta_d() << '\n';
   ostr << "zeta_l = " << get_zeta_l() << '\n';

   ostr << "yuh =\n"  << get_yuh()  << '\n';
   ostr << "yuH =\n"  << get_yuH()  << '\n';
   ostr << "yuA =\n"  << get_yuA()  << '\n';
   ostr << "yuHp =\n" << get_yuHp() << '\n';
   ostr << "ydh =\n"  << get_ydh()  << '\n';
   ostr << "ydH =\n"  << get_ydH()  << '\n';
   ostr << "ydA =\n"  << get_ydA()  << '\n';
   ostr << "ydHp =\n" << get_ydHp() << '\n';
   ostr << "ylh =\n"  << get_ylh()  << '\n';
   ostr << "ylH =\n"  << get_ylH()  << '\n';
   ostr << "ylA =\n"  << get_ylA()  << '\n';
   ostr << "ylHp =\n" << get_ylHp() << '\n';

   ostr << get_sm();
}

const char* THDM::yukawa_type_to_string() const
{
   switch (yukawa_type) {
      case thdm::Yukawa_type::type_1:
         return "Type I";
         break;
      case thdm::Yukawa_type::type_2:
         return "Type II";
         break;
      case thdm::Yukawa_type::type_X:
         return "Type X";
         break;
      case thdm::Yukawa_type::type_Y:
         return "Type Y";
         break;
      case thdm::Yukawa_type::aligned:
         return "Aligned";
         break;
      case thdm::Yukawa_type::general:
         return "General";
         break;
   }
   throw ESetupError("Bug: unhandled case in yukawa_type_to_string.");
}

void THDM::validate() const
{
   // check zeta_f
   if (yukawa_type == thdm::Yukawa_type::type_1 ||
       yukawa_type == thdm::Yukawa_type::type_2 ||
       yukawa_type == thdm::Yukawa_type::type_X ||
       yukawa_type == thdm::Yukawa_type::type_Y) {
      if (zeta_u != 0) {
         WARNING("Value of zeta_u = " << zeta_u << " ignored, because Yukawa type is " << yukawa_type_to_string());
      }
      if (zeta_d != 0) {
         WARNING("Value of zeta_d = " << zeta_d << " ignored, because Yukawa type is " << yukawa_type_to_string());
      }
      if (zeta_l != 0) {
         WARNING("Value of zeta_l = " << zeta_l << " ignored, because Yukawa type is " << yukawa_type_to_string());
      }
   }

   // check X_f
   if (yukawa_type == thdm::Yukawa_type::type_1 ||
       yukawa_type == thdm::Yukawa_type::type_2 ||
       yukawa_type == thdm::Yukawa_type::type_X ||
       yukawa_type == thdm::Yukawa_type::type_Y ||
       yukawa_type == thdm::Yukawa_type::aligned) {
      if (get_Xu().cwiseAbs().maxCoeff() != 0) {
         WARNING("Value of Xu ignored, because Yukawa type is " << yukawa_type_to_string());
      }
      if (get_Xd().cwiseAbs().maxCoeff() != 0) {
         WARNING("Value of Xd ignored, because Yukawa type is " << yukawa_type_to_string());
      }
      if (get_Xl().cwiseAbs().maxCoeff() != 0) {
         WARNING("Value of Xl ignored, because Yukawa type is " << yukawa_type_to_string());
      }
   }
}

std::ostream& operator<<(std::ostream& ostr, const THDM& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace gm2calc
