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

namespace gm2calc {

class MSSMNoFV_onshell;

/// calculates full 1-loop SUSY contributions to (g-2) in the MSSM (w/ tan(beta) resummation)
double calculate_amu_1loop(const MSSMNoFV_onshell&);

/// calculates full 1-loop SUSY contributions to (g-2) in the MSSM (no tan(beta) resummation)
double calculate_amu_1loop_non_tan_beta_resummed(const MSSMNoFV_onshell&);

// === routines for individual 1-loop contributions ===

/// 1-loop neutralino contribution
double amuChi0(const MSSMNoFV_onshell&);

/// 1-loop chargino contribution
double amuChipm(const MSSMNoFV_onshell&);

// === approximations ===

// main routines

/// 1-loop leading log approximation
double amu1Lapprox(const MSSMNoFV_onshell&);

/// 1-loop leading log approximation w/o explicit tan(beta) resummation
double amu1Lapprox_non_tan_beta_resummed(const MSSMNoFV_onshell&);

// routines for individual contributions

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

// === resummation ===

double tan_beta_cor(const MSSMNoFV_onshell&);

} // namespace gm2calc

#endif
