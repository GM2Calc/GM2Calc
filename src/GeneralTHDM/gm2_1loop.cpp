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

#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/GeneralTHDM.hpp"
#include "GeneralTHDM/gm2_1loop_helpers.hpp"
#include <cmath>

/**
 * \file gm2_1loop.cpp
 *
 * Contains functions necessary to calculate the THDM
 * contributions for g-2 at the 1-loop level.
 */

namespace gm2calc {

namespace {

const double sqrt2 = 1.4142135623730950; // Sqrt[2]
double sqr(double x) noexcept { return x*x; }

} // anonymous namespace

/**
 * Calculates full 1-loop contribution to a_mu in the general THDM.
 *
 * @todo(alex) check convention for ylHp
 *
 * @param model THDM model parameters, masses and mixings
 * @return 1-loop contribution to a_mu
 */
double calculate_amu_1loop(const GeneralTHDM& model)
{
   const double alpha_h = model.get_alpha_h();
   const double beta = model.get_beta();
   const double tb = model.get_tan_beta();
   const double cb = 1./std::sqrt(1.0 + sqr(tb));
   const double sba = std::sin(beta - alpha_h);
   const double cba = std::cos(beta - alpha_h);
   const double v = model.get_v();
   const double mm = model.get_MFe(1);
   const Eigen::Matrix<double,3,3> id = Eigen::Matrix<double,3,3>::Identity();

   const Eigen::Matrix<double,3,3> zetal =
      model.get_Xe()*v/(cb*sqrt2*mm) + model.get_zeta_l()*id;

   general_thdm::THDM_1L_parameters pars;
   pars.alpha_em = model.get_alpha_em();
   pars.mm = model.get_MFe(1);
   pars.mw = model.get_MVWm();
   pars.mz = model.get_MVZ();
   pars.mhSM = model.get_physical().MhSM;
   pars.mA = model.get_MAh(1);
   pars.mHp = model.get_MHm(1);
   pars.ml = model.get_MFe();
   pars.mv = model.get_MFv();
   pars.mh = model.get_Mhh();
   pars.ylh = sba*id + cba*zetal; // Eq.(18), arxiv:1607.06292
   pars.ylH = cba*id - sba*zetal; // Eq.(18), arxiv:1607.06292
   pars.ylA = -zetal;             // Eq.(18), arxiv:1607.06292
   pars.ylHp = pars.ylA; // @todo(alex) check convention

   return general_thdm::amu1L(pars);
}

} // namespace gm2calc
