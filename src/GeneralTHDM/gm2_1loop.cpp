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
#include <cmath>

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

/// Eq.(28), arxiv:1607.06292
double Fh(double x) noexcept
{
   return -std::log(x) - 7.0/6;
}

/// Eq.(29), arxiv:1607.06292
double FA(double x) noexcept
{
   return std::log(x) + 11.0/6;
}

/// Eq.(30), arxiv:1607.06292
double FHp(double x) noexcept
{
   return -1.0/6;
}

} // anonymous namespace

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
      + sqr(pars.ylHp(1,1)) * mm2/mHp2 * FHp(mm2/mHp2) // @todo(check) Yukawa prefactor
      // subtract SM contribution
      - mm2/mhSM2 * Fh(mm2/mhSM2);

   // Eq.(33), arxiv:1607.06292
   const auto GF = pi*pars.alpha/(sqrt2*sw2*mw2)*(1 + delta_r);
   const auto pref = GF*mm2/(4*sqrt2*pi2);

   return pref*res;
}

double amu1L(const THDM_1L_parameters&) noexcept
{
   return 0.0;
}

} // namespace general_thdm

} // namespace gm2calc
