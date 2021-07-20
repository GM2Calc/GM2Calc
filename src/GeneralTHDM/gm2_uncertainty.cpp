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

#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/GeneralTHDM.hpp"
#include "gm2_1loop_helpers.hpp"

#include <cmath>

/**
 * \file gm2_uncertainty.cpp
 *
 * Contains functions necessary to calculate theory uncertainty
 * associated with amu.
 */

namespace gm2calc {

/**
 * Calculates uncertainty associated with amu(0-loop)
 *
 * The estimated uncertainty is the magnitude amu(1-loop) plus
 * amu(2-loop).
 *
 * @param model model parameters
 *
 * @return uncertainty for amu(0-loop)
 */
double calculate_uncertainty_amu_0loop(const GeneralTHDM& model)
{
   const double amu_1L = calculate_amu_1loop(model);
   const double amu_2L = calculate_amu_2loop(model);

   return std::abs(amu_1L) + std::abs(amu_2L);
}

/**
 * Calculates uncertainty associated with amu(1-loop)
 *
 * The estimated uncertainty is the sum of magnitude amu(2-loop) and
 * the 2-loop uncertainty, calculated by
 * calculate_uncertainty_amu_2loop().
 *
 * @param model model parameters
 *
 * @return uncertainty for amu(1-loop)
 */
double calculate_uncertainty_amu_1loop(const GeneralTHDM& model)
{
   const double amu_2L = calculate_amu_2loop(model);
   const double delta_amu_2L = calculate_uncertainty_amu_2loop(model);

   return std::abs(amu_2L) + std::abs(delta_amu_2L);
}

/**
 * Calculates uncertainty associated with amu(2-loop)
 *
 * Takes into account the neglected two-loop contribution
 * \f$a_\mu^{\Delta r\text{-shift}}\f$, Eq.(34) arxiv:1607.06292.
 *
 * @param model model parameters
 *
 * @return uncertainty for amu(2-loop)
 */
double calculate_uncertainty_amu_2loop(const GeneralTHDM& model)
{
   const double delta_r_shift_max = 2e-12;

   const double alpha_em = model.get_alpha_em();
   const double mHp = model.get_MHm(1);
   const double mz = model.get_MVZ();
   const double delta_alpha_em_shift
      = calculate_amu_2loop(model)*general_thdm::delta_alpha(alpha_em, mHp, mz);

   return delta_r_shift_max + std::abs(delta_alpha_em_shift);
}

} // namespace gm2calc
