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

#ifndef GM2_1LOOP_H
#define GM2_1LOOP_H

#include <cmath>
#include <Eigen/Core>

namespace gm2calc {

class MSSMNoFV_onshell;

static const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);

/// calculates full 1-loop SUSY contributions to (g-2) in the MSSM (w/ tan(beta) resummation)
double calculate_amu_1loop(const MSSMNoFV_onshell&);

/// calculates full 1-loop SUSY contributions to (g-2) in the MSSM (no tan(beta) resummation)
double calculate_amu_1loop_non_tan_beta_resummed(const MSSMNoFV_onshell&);

/// 1-loop neutralino contribution
double amuChi0(const MSSMNoFV_onshell&);

/// 1-loop chargino contribution
double amuChipm(const MSSMNoFV_onshell&);

// === approximations ===

/// 1-loop leading log approximation
double amu1Lapprox(const MSSMNoFV_onshell&);

/// 1-loop wino--Higgsino, muon-sneutrino leading log approximation
double amuWHnu(const MSSMNoFV_onshell&);
/// 1-loop wino--Higgsino, left-handed smuon leading log approximation
double amuWHmuL(const MSSMNoFV_onshell&);
/// 1-loop bino--Higgsino, left-handed smuon leading log approximation
double amuBHmuL(const MSSMNoFV_onshell&);
/// 1-loop bino--Higgsino, right-handed smuon leading log approximation
double amuBHmuR(const MSSMNoFV_onshell&);
/// 1-loop bino, left-handed smuon--right-handed smuon leading log approximation
double amuBmuLmuR(const MSSMNoFV_onshell&);

// === resummations ===

double tan_beta_cor(const MSSMNoFV_onshell&);
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
