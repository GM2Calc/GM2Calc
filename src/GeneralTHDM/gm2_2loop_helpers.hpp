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

#ifndef GM2_GENERALTHDM_2LOOP_HELPERS_HPP
#define GM2_GENERALTHDM_2LOOP_HELPERS_HPP

#include <Eigen/Core>

namespace gm2calc {

// === 2-loop bosonic contributions ===

// routines for sub-expressions

double amu2L_B_EWadd(double eta, double zetal);
double amu2L_B_nonYuk();
double amu2L_B_Yuk();
double amu2L_F();

} // namespace gm2calc

#endif
