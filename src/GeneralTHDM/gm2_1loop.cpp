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

} // anonymous namespace

double amu1L(const THDM_1L_parameters& pars) noexcept
{
   const auto mm2 = pars.mm*pars.mm;
   const auto mw2 = pars.mw*pars.mw;
   const auto mz2 = pars.mz*pars.mz;
   const auto cw2 = mw2/mz2;
   const auto sw2 = 1 - cw2;
   const auto delta_r = 0.0; // @todo(alex)

   // Eq.(27), arxiv:1607.06292
   const auto res = 0.0; // @todo(alex)

   // Eq.(33), arxiv:1607.06292
   const auto GF = pi*pars.alpha/sqrt2*sw2*mw2*(1 + delta_r);
   const auto pref = GF*mm2/(4*sqrt2*pi2);

   return pref*res;
}

} // namespace general_thdm

} // namespace gm2calc
