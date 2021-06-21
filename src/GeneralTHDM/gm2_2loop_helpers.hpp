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

/// parameters to be passed to the fermionic contribution functions
struct THDM_F_parameters {
   double alpha{}; ///< alpha_em
   double mw{};    ///< W boson mass
   double mz{};    ///< Z boson mass
   Eigen::Matrix<double,2,1> mh{Eigen::Matrix<double,2,1>::Zero()}; ///< CP-even Higgs bosons mass
   double mhSM{};  ///< SM Higgs boson mass
   double mA{};    ///< CP-odd Higgs boson mass
   double mHp{};   ///< charged Higgs boson mass
   Eigen::Matrix<double,3,1> ml{Eigen::Matrix<double,3,1>::Zero()}; ///< down-type lepton masses
   Eigen::Matrix<double,3,1> mu{Eigen::Matrix<double,3,1>::Zero()}; ///< up-type quark masses
   Eigen::Matrix<double,3,1> md{Eigen::Matrix<double,3,1>::Zero()}; ///< down-type quark masses
   double yuh{};   ///< Y_f^h coefficients with f=u and S=h
   double ydh{};   ///< Y_f^h coefficients with f=d and S=h
   double ylh{};   ///< Y_f^h coefficients with f=l and S=h
   double yuH{};   ///< Y_f^H coefficients with f=u and S=H
   double ydH{};   ///< Y_f^H coefficients with f=d and S=H
   double ylH{};   ///< Y_f^H coefficients with f=l and S=H
   double yuA{};   ///< Y_f^A coefficients with f=u and S=A
   double ydA{};   ///< Y_f^A coefficients with f=d and S=A
   double ylA{};   ///< Y_f^A coefficients with f=l and S=A
};

// === 2-loop bosonic contributions ===

// routines for sub-expressions

double amu2L_B_EWadd(double eta, double zetal);
double amu2L_B_nonYuk();
double amu2L_B_Yuk();

// === 2-loop fermionic contributions ===

double amu2L_F(const THDM_F_parameters&);

} // namespace gm2calc

#endif
