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

#include "GeneralTHDM/gm2_2loop_helpers.hpp"

/**
 * \file gm2_2loop.cpp
 *
 * Contains functions necessary to calculate the THDM contributions
 * for g-2 at the 2-loop level.
 */

namespace gm2calc {

/**
 * \fn amu2L_B_EWadd
 *
 * Calculates 2-loop bosonic pure electroweak contributions.
 *
 * Eq (49), arxiv:1607:06292
 */
double amu2L_B_EWadd(double eta, double zetal)
{
   return 2.3e-11 * eta * zetal;
}

double amu2L_B_nonYuk()
{
   return 0.0;
}

double amu2L_B_Yuk()
{
   return 0.0;
}

} // namespace gm2calc