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

/**
 * \file gm2_uncertainty.cpp
 *
 * Contains functions necessary to calculate theory uncertainty
 * associated with amu.
 */

namespace gm2calc {

/**
 * Calculates uncertainty associated with amu(2-loop) using Eq (84),
 * arXiv:hep-ph/0609168 .
 *
 * @note Eq (84), arXiv:hep-ph/0609168 takes into account photonic
 * 2-loop corrections and 2-loop fermion/sfermion corrections, which
 * were unknown at the time arXiv:hep-ph/0609168 was published.
 * However, GM2Calc includes the complete photonic 2-loop corrections
 * and approximations for the 2-loop fermion/sfermion corrections.
 * For this reason it is debatable whether Eq. (84) (especially the
 * finite offset) is applicable for the 2-loop calculation of GM2Calc.
 *
 * @param amu_1L amu at 1-loop level
 *
 * @return uncertainty for amu(2-loop)
 */
double calculate_uncertainty_amu_2loop(double amu_1L)
{
   return 0.02 * amu_1L + 2.5e-10;
}

} // namespace gm2calc
