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

#include "GeneralTHDM/gm2_1loop_helpers.hpp"
#include "gm2_ffunctions.hpp"
#include <cmath>
#include <complex>

/**
 * \file gm2_1loop.cpp
 *
 * Contains functions necessary to calculate the THDM
 * contributions for g-2 at the 1-loop level.
 */

namespace gm2calc {

namespace general_thdm {

namespace {

const double pi = 3.1415926535897932;
const double pi2 = 9.8696044010893586; // Pi^2
const double sqrt2 = 1.4142135623730950; // Sqrt[2]

double sqr(double x) noexcept { return x*x; }

double conj(double x) noexcept { return x; }

/// Eq.(28), arxiv:1607.06292
double Fh(double x) noexcept
{
   return F1C(x)/12 + F2C(x)/3;
}

/// Eq.(29), arxiv:1607.06292
double FA(double x) noexcept
{
   return F1C(x)/12 - F2C(x)/3;
}

/// Eq.(30), arxiv:1607.06292
double FHp(double x) noexcept
{
   return -F1N(x)/12;
}

double AS(int gen, const Eigen::Matrix<double,3,1>& ml, double mS2, const Eigen::Matrix<double,3,3>& y) noexcept
{
   const auto x = sqr(ml(gen))/mS2;
   const double y2 = conj(y(gen, 1))*conj(y(1, gen));

   return
      + 1.0/12*F1C(x)*(std::norm(y(gen, 1)) + std::norm(y(1, gen)))
      + 2.0/3*ml(gen)/ml(1)*F2C(x)*y2;
}

double AA(int gen, const Eigen::Matrix<double,3,1>& ml, double mS2, const Eigen::Matrix<double,3,3>& y) noexcept
{
   const auto x = sqr(ml(gen))/mS2;
   const double y2 = conj(y(gen, 1))*conj(y(1, gen));

   return
      + 1.0/12*F1C(x)*(std::norm(y(gen, 1)) + std::norm(y(1, gen)))
      - 2.0/3*ml(gen)/ml(1)*F2C(x)*y2;
}

double AHp(int gen, const Eigen::Matrix<double,3,1>& mv, double mS2, const Eigen::Matrix<double,3,3>& y) noexcept
{
   return std::norm(y(gen, 1))/24*(
      F1N(sqr(mv(1))/mS2) + F1N(sqr(mv(gen))/mS2));
}


/// 1-loop THDM contribution to \f$\Delta r\f$
double delta_alpha(double alpha, double mHp, double q) noexcept
{
   return -alpha/(6*pi)*std::log(std::abs(mHp/q));
}

/// 1-loop THDM contribution to \f$\Delta\rho\f$
double delta_rho() noexcept
{
   return 0.0; // @todo(alex)
}

/// 1-loop THDM contribution to \f$\Delta r_{\text{rem}}\f$
double delta_rem() noexcept
{
   return 0.0; // @todo(alex)
}

} // anonymous namespace


/**
 * Calculation of the 1-loop THDM contribution to \f$\Delta r\f$
 * from arxiv:1211.0311
 */
double delta_r(const THDM_delta_r_parameters& pars) noexcept
{
   const auto mw2 = sqr(pars.mw);
   const auto mz2 = sqr(pars.mz);
   const auto cw2 = mw2/mz2;
   const auto sw2 = 1 - cw2;
   const auto d_alpha = delta_alpha(pars.alpha, pars.mHp, pars.mw);
   const auto d_rho = delta_rho();
   const auto d_rem = delta_rem();

   // Eq.(11), arxiv:1211.0311
   return d_alpha - cw2/sw2*d_rho + d_rem;
}

/**
 * Approximation for 1-loop contribution
 * Eq (27) from arxiv:1607.06292
 */
double amu1L_approx(const THDM_1L_parameters& pars) noexcept
{
   const auto mm2 = sqr(pars.mm);
   const auto mw2 = sqr(pars.mw);
   const auto mz2 = sqr(pars.mz);
   const auto cw2 = mw2/mz2;
   const auto sw2 = 1 - cw2;
   const auto mhSM2 = sqr(pars.mhSM);
   const auto mh2 = sqr(pars.mh(0));
   const auto mH2 = sqr(pars.mh(1));
   const auto mA2 = sqr(pars.mA);
   const auto mHp2 = sqr(pars.mHp);
   const auto delta_r = 0.0; // @todo(alex)

   // Eq.(27), arxiv:1607.06292
   const auto res =
      + sqr(pars.ylh(1,1)) * mm2/mh2 * Fh(mm2/mh2)
      + sqr(pars.ylH(1,1)) * mm2/mH2 * Fh(mm2/mH2)
      + sqr(pars.ylA(1,1)) * mm2/mA2 * FA(mm2/mA2)
      + sqr(pars.ylHp(1,1)) * mm2/mHp2 * FHp(mm2/mHp2)
      // subtract SM contribution
      - mm2/mhSM2 * Fh(mm2/mhSM2);

   // Eq.(33), arxiv:1607.06292
   const auto GF = pi*pars.alpha/(sqrt2*sw2*mw2)*(1 + delta_r);
   const auto pref = GF*mm2/(4*sqrt2*pi2);

   return pref*res;
}

/**
 * 1-loop contribution
 *
 * @note The CP-odd couplings are assumed to be real
 * @todo(alex) check prefactor 1/2 for neutral contributions
 * @todo(alex) check prefactor (-1) for charged contributions
 */
double amu1L(const THDM_1L_parameters& pars) noexcept
{
   const auto mm2 = sqr(pars.mm);
   const auto mw2 = sqr(pars.mw);
   const auto mz2 = sqr(pars.mz);
   const auto cw2 = mw2/mz2;
   const auto sw2 = 1 - cw2;
   const auto mhSM2 = sqr(pars.mhSM);
   const auto mh2 = sqr(pars.mh(0));
   const auto mH2 = sqr(pars.mh(1));
   const auto mA2 = sqr(pars.mA);
   const auto mHp2 = sqr(pars.mHp);
   const auto delta_r = 0.0; // @todo(alex)
   Eigen::Matrix<double,3,3> ylhSM{Eigen::Matrix<double,3,3>::Zero()};
   ylhSM(1,1) = 1.0;

   double res = 0.0;

   for (int g = 0; g < 3; ++g) {
      res += 0.5 * mm2/mh2 * AS(g, pars.ml, mh2, pars.ylh);
      res += 0.5 * mm2/mH2 * AS(g, pars.ml, mH2, pars.ylH);
      res += 0.5 * mm2/mA2 * AA(g, pars.ml, mA2, pars.ylA);
      res += -mm2/mHp2 * AHp(g, pars.mv, mHp2, pars.ylHp);
      // subtract SM contribution
      res -= 0.5 * mm2/mhSM2 * AS(g, pars.ml, mhSM2, ylhSM);
   }

   const auto GF = pi*pars.alpha/(sqrt2*sw2*mw2)*(1 + delta_r);
   const auto pref = GF*mm2/(4*sqrt2*pi2);

   return pref*res;
}

} // namespace general_thdm

} // namespace gm2calc
