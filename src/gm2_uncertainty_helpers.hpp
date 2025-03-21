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

#ifndef GM2_UNCERTAINTY_HELPERS_HPP
#define GM2_UNCERTAINTY_HELPERS_HPP

namespace gm2calc {

class MSSMNoFV_onshell;
class THDM;

/// calculates uncertainty for amu(0-loop)
double calculate_uncertainty_amu_0loop(const THDM& /* model */, double /* amu_1L */, double /* amu_2L */);

/// calculates uncertainty for amu(1-loop)
double calculate_uncertainty_amu_1loop(const THDM& /* model */, double /* amu_1L */, double /* amu_2L */);

/// calculates uncertainty for amu(2-loop)
double calculate_uncertainty_amu_2loop(const THDM& /* model */, double /* amu_1L */, double /* amu_2L */);

/// calculates uncertainty for amu(0-loop) w/ tan(beta) resummation
double calculate_uncertainty_amu_0loop(const MSSMNoFV_onshell& /* model */, double /* amu_1L*/);

/// calculates uncertainty for amu(1-loop) w/ tan(beta) resummation
double calculate_uncertainty_amu_1loop(const MSSMNoFV_onshell& /* model */, double /* amu_2L */);

} // namespace gm2calc

#endif
