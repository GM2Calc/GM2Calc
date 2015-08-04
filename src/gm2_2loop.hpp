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

#ifndef GM2_2LOOP_H
#define GM2_2LOOP_H

#include <Eigen/Core>

namespace flexiblesusy {
namespace gm2calc {

class MSSMNoFV_onshell;

// === 2-loop fermion/sfermion approximations ===

double amu2LFSfapprox(const MSSMNoFV_onshell&);
double amu2LFSfapprox_non_tan_beta_resummed(const MSSMNoFV_onshell&);

double amuWHnu2L(const MSSMNoFV_onshell&);
double amuWHmuL2L(const MSSMNoFV_onshell&);
double amuBHmuL2L(const MSSMNoFV_onshell&);
double amuBHmuR2L(const MSSMNoFV_onshell&);
double amuBmuLmuR2L(const MSSMNoFV_onshell&);

double LogNorm(const MSSMNoFV_onshell&);

double Deltag1(const MSSMNoFV_onshell&);
double DeltaYukHiggsino(const MSSMNoFV_onshell&);
double DeltaYukBinoHiggsino(const MSSMNoFV_onshell&);
double Deltag2(const MSSMNoFV_onshell&);
double DeltaYukWinoHiggsino(const MSSMNoFV_onshell&);
double DeltaTanBeta(const MSSMNoFV_onshell&);


// === photonic 2-loop corrections ===

double amuChipmPhotonic(const MSSMNoFV_onshell&);
double amuChi0Photonic(const MSSMNoFV_onshell&);

// === SUSY 2L(a) diagrams ===

double amua2LSferm(const MSSMNoFV_onshell&);
double amua2LCha(const MSSMNoFV_onshell&);

double tan_alpha(const MSSMNoFV_onshell&);
Eigen::Matrix<std::complex<double>,3,3> lambda_mu_cha(const MSSMNoFV_onshell&);
Eigen::Matrix<std::complex<double>,2,2> lambda_stop(const MSSMNoFV_onshell&);
Eigen::Matrix<std::complex<double>,2,2> lambda_sbot(const MSSMNoFV_onshell&);
Eigen::Matrix<std::complex<double>,2,2> lambda_stau(const MSSMNoFV_onshell&);

} // namespace gm2calc
} // namespace flexiblesusy

#endif
