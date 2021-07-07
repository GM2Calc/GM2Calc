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

#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/GeneralTHDM.hpp"
#include "GeneralTHDM/gm2_2loop_helpers.hpp"

namespace gm2calc {

/**
 * Calculates full 2-loop contribution to a_mu in the general THDM.
 *
 * @todo(alex) to be implemented
 *
 * @param model THDM model parameters, masses and mixings
 * @return 2-loop contribution to a_mu
 */
double calculate_amu_2loop(const GeneralTHDM& model)
{
   general_thdm::THDM_B_parameters pars_b;
   pars_b.alpha = model.get_alpha_em();
   pars_b.mm = model.get_MFe(1);
   pars_b.mw = model.get_MVWm();
   pars_b.mz = model.get_MVZ();
   pars_b.mhSM = model.get_physical().MhSM;
   pars_b.mA = model.get_MAh(1);
   pars_b.mHp = model.get_MHm(1);
   pars_b.mh = model.get_Mhh();
   pars_b.tb = model.get_tan_beta();
   // pars_b.zetal = ?;
   pars_b.eta = model.get_eta();
   pars_b.lambda5 = model.get_LambdaFive();

   general_thdm::THDM_F_parameters pars_f;
   pars_f.alpha = model.get_alpha_em();
   pars_f.mm = model.get_MFe(1);
   pars_f.mw = model.get_MVWm();
   pars_f.mz = model.get_MVZ();
   pars_f.mhSM = model.get_physical().MhSM;
   pars_f.mA = model.get_MAh(1);
   pars_f.mHp = model.get_MHm(1);
   pars_f.mh = model.get_Mhh();
   pars_f.ml = model.get_MFe();
   pars_f.mu = model.get_MFu();
   pars_f.md = model.get_MFd();
   // pars_b.yuS = ?;
   // pars_b.ydS = ?;
   // pars_b.ylS = ?;

   return amu2L_B(pars_b) + amu2L_F(pars_f);
}

} // namespace gm2calc
