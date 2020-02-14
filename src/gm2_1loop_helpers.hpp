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

#ifndef GM2_1LOOP_HELPERS_H
#define GM2_1LOOP_HELPERS_H

#include <Eigen/Core>

namespace gm2calc {

class MSSMNoFV_onshell;

// === resummations ===

double delta_mu_correction(const MSSMNoFV_onshell&);
double delta_tau_correction(const MSSMNoFV_onshell&);
double delta_bottom_correction(const MSSMNoFV_onshell&);

// === couplings ===

Eigen::Array<double,2,1> AAC(const MSSMNoFV_onshell&);
Eigen::Matrix<double,4,2> AAN(const MSSMNoFV_onshell&);
Eigen::Array<double,2,1> BBC(const MSSMNoFV_onshell&);
Eigen::Matrix<double,4,2> BBN(const MSSMNoFV_onshell&);

/// squared neutralino smuon mass ratio
Eigen::Matrix<double,4,2> x_im(const MSSMNoFV_onshell&);
/// squared chargino muon-sneutrino mass ratio
Eigen::Array<double,2,1> x_k(const MSSMNoFV_onshell&);

} // namespace gm2calc

#endif
