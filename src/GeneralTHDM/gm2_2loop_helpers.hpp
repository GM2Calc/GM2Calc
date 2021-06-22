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

namespace general_thdm {

/// parameters to be passed to the bosonic contribution functions
struct THDM_B_parameters {
   double alpha{};   ///< alpha_em
   double mw{};      ///< W boson mass
   double mz{};      ///< Z boson mass
   double mm{};      ///< muon mass for prefactor
   double mhSM{};    ///< SM Higgs boson mass
   double mA{};      ///< CP-odd Higgs boson mass
   double mHp{};     ///< charged Higgs boson mass
   Eigen::Matrix<double,2,1> mh{Eigen::Matrix<double,2,1>::Zero()}; ///< CP-even Higgs bosons mass
   double tb{};      ///< tan(beta)
   double zetal{};   ///< zeta_l
   double eta{};     ///< eta
   double lambda5{}; ///< Lambda_5
};

/// parameters to be passed to the fermionic contribution functions
struct THDM_F_parameters {
   double alpha{}; ///< alpha_em
   double mw{};    ///< W boson mass
   double mz{};    ///< Z boson mass
   double mhSM{};  ///< SM Higgs boson mass
   double mA{};    ///< CP-odd Higgs boson mass
   double mHp{};   ///< charged Higgs boson mass
   Eigen::Matrix<double,2,1> mh{Eigen::Matrix<double,2,1>::Zero()}; ///< CP-even Higgs bosons mass
   Eigen::Matrix<double,3,1> ml{Eigen::Matrix<double,3,1>::Zero()}; ///< down-type lepton masses
   Eigen::Matrix<double,3,1> mu{Eigen::Matrix<double,3,1>::Zero()}; ///< up-type quark masses
   Eigen::Matrix<double,3,1> md{Eigen::Matrix<double,3,1>::Zero()}; ///< down-type quark masses
   Eigen::Matrix<double,3,3> yuS{Eigen::Matrix<double,3,3>::Zero()};///< Y_f^S coefficients with f={u,c,t} and S={h,H,A}
   Eigen::Matrix<double,3,3> ydS{Eigen::Matrix<double,3,3>::Zero()};///< Y_f^S coefficients with f={d,s,b} and S={h,H,A}
   Eigen::Matrix<double,3,3> ylS{Eigen::Matrix<double,3,3>::Zero()};///< Y_f^S coefficients with f={e,m,τ} and S={h,H,A}
};

// === 2-loop bosonic contributions ===

double amu2L_B(const THDM_B_parameters&);

// routines for sub-expressions

double amu2L_B_EWadd(const THDM_B_parameters&);
double amu2L_B_nonYuk(const THDM_B_parameters&);
double amu2L_B_Yuk(const THDM_B_parameters&);

// === 2-loop fermionic contributions ===

double amu2L_F(const THDM_F_parameters&);

// routines for sub-expressions

double amu2L_F_charged(const THDM_F_parameters&);
double amu2L_F_neutral(const THDM_F_parameters&);

} // namespace general_thdm

} // namespace gm2calc

#endif
