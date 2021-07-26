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

#include "gm2calc/SM.h"
#include "gm2calc/SM.hpp"

extern "C"
{

/**
 * @brief Allocate a new SM.
 * @return pointer to SM
 */
gm2calc_SM* gm2calc_generalthdm_sm_new()
{
   return reinterpret_cast<gm2calc_SM*>(new gm2calc::SM());
}

/**
 * @brief Deletes a SM.
 * @return pointer to SM
 */
void gm2calc_generalthdm_sm_free(gm2calc_SM* sm)
{
   delete reinterpret_cast<gm2calc::SM*>(sm);
}

void gm2calc_generalthdm_sm_set_alpha_em_0(gm2calc_SM* sm, double val)
{
   reinterpret_cast<gm2calc::SM*>(sm)->set_alpha_em_0(val);
}

void gm2calc_generalthdm_sm_set_alpha_em_mz(gm2calc_SM* sm, double val)
{
   reinterpret_cast<gm2calc::SM*>(sm)->set_alpha_em_mz(val);
}

void gm2calc_generalthdm_sm_set_mh(gm2calc_SM* sm, double val)
{
   reinterpret_cast<gm2calc::SM*>(sm)->set_mh(val);
}

void gm2calc_generalthdm_sm_set_mw(gm2calc_SM* sm, double val)
{
   reinterpret_cast<gm2calc::SM*>(sm)->set_mw(val);
}

void gm2calc_generalthdm_sm_set_mz(gm2calc_SM* sm, double val)
{
   reinterpret_cast<gm2calc::SM*>(sm)->set_mz(val);
}

void gm2calc_generalthdm_sm_set_mu(gm2calc_SM* sm, unsigned i, double val)
{
   reinterpret_cast<gm2calc::SM*>(sm)->set_mu(i, val);
}

void gm2calc_generalthdm_sm_set_md(gm2calc_SM* sm, unsigned i, double val)
{
   reinterpret_cast<gm2calc::SM*>(sm)->set_md(i, val);
}

void gm2calc_generalthdm_sm_set_ml(gm2calc_SM* sm, unsigned i, double val)
{
   reinterpret_cast<gm2calc::SM*>(sm)->set_ml(i, val);
}

void gm2calc_generalthdm_sm_set_ckm_real(gm2calc_SM* sm, unsigned i, unsigned j, double re)
{
   auto m = reinterpret_cast<gm2calc::SM*>(sm);
   const auto im = std::imag(m->get_ckm(i, j));
   m->set_ckm(i, j, std::complex<double>(re, im));
}

void gm2calc_generalthdm_sm_set_ckm_imag(gm2calc_SM* sm, unsigned i, unsigned j, double im)
{
   auto m = reinterpret_cast<gm2calc::SM*>(sm);
   const auto re = std::imag(m->get_ckm(i, j));
   m->set_ckm(i, j, std::complex<double>(re, im));
}

void gm2calc_generalthdm_sm_set_ckm_from_wolfenstein(gm2calc_SM* sm, double lambdaW, double aCkm, double rhobar, double etabar)
{
   reinterpret_cast<gm2calc::SM*>(sm)->set_ckm_from_wolfenstein(lambdaW, aCkm, rhobar, etabar);
}

void gm2calc_generalthdm_sm_set_ckm_from_angles(gm2calc_SM* sm, double theta_12, double theta_13, double theta_23, double delta)
{
   reinterpret_cast<gm2calc::SM*>(sm)->set_ckm_from_angles(theta_12, theta_13, theta_23, delta);
}

} // extern "C"
