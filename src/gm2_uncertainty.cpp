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

#include "gm2_uncertainty.hpp"
#include "gm2_2loop.hpp"
#include <cmath>

/**
 * \file gm2_uncertainty.cpp
 *
 * Contains functions necessary to calculate theory uncertainty
 * associated with amu.
 */

namespace gm2calc {

/**
 * Calculates uncertainty associated with amu(2-loop) using Eq (4).
 *
 * Eq. (4) takes into account the unknown two-loop contributions and
 * the employed approximation for the 2L(a) contributions.
 *
 * @param model model parameters
 *
 * @return uncertainty for amu(2-loop)
 */
double calculate_uncertainty_amu_2loop(const MSSMNoFV_onshell& model)
{
   using std::abs;
   const double amu_2La_Cha = amu2LaCha(model);
   const double amu_2La_Sferm = amu2LaSferm(model);

   return 2.3e-10 + 0.3 * (abs(amu_2La_Cha) + abs(amu_2La_Sferm));
}

} // namespace gm2calc
