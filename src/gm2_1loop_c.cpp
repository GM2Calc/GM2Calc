/* ====================================================================
 * This file is part of GM2Calc.
 *
 * GM2Calc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * GM2Calc is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GM2Calc.  If not, see
 * <http://www.gnu.org/licenses/>.
 * ==================================================================== */

#include "gm2_1loop.h"
#include "gm2_1loop.hpp"
#include "MSSMNoFV_onshell.hpp"

extern "C" {

/** calculates full 1-loop SUSY contributions to (g-2) in the MSSM (w/ tan(beta) resummation) */
double gm2calc_calculate_amu_1loop(MSSMNoFV_onshell* model)
{
   return gm2calc::calculate_amu_1loop(
      *reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model));
}

/** calculates full 1-loop SUSY contributions to (g-2) in the MSSM (no tan(beta) resummation) */
double gm2calc_calculate_amu_1loop_non_tan_beta_resummed(MSSMNoFV_onshell* model)
{
   return gm2calc::calculate_amu_1loop_non_tan_beta_resummed(
      *reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model));
}

/* === routines for individual 1-loop contributions === */

/** 1-loop neutralino contribution */
double gm2calc_amuChi0(MSSMNoFV_onshell* model)
{
   return gm2calc::amuChi0(
      *reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model));
}

/** 1-loop chargino contribution */
double gm2calc_amuChipm(MSSMNoFV_onshell* model)
{
   return gm2calc::amuChipm(
      *reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model));
}

} /* extern "C" */
