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

#include "gm2_2loop.h"
#include "gm2_2loop.hpp"

/**
 * @file gm2_2loop_c.cpp
 * @brief contains definitions of C interface functions for 2-loop calculation
 *
 * This file contains the definitions for the C interface functions
 * used to calculate \f$a_\mu\f$ at the 2-loop level.
 */

extern "C" {

/** calculates best 2-loop SUSY contributions to a_mu in the MSSM (with tan(beta) resummation) */
double gm2calc_calculate_amu_2loop(const MSSMNoFV_onshell* model)
{
   return gm2calc::calculate_amu_2loop(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

/** calculates best 2-loop SUSY contributions to a_mu in the MSSM (no tan(beta) resummation) */
double gm2calc_calculate_amu_2loop_non_tan_beta_resummed(const MSSMNoFV_onshell* model)
{
   return gm2calc::calculate_amu_2loop_non_tan_beta_resummed(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

/* === 2-loop fermion/sfermion approximations === */

double gm2calc_amu2LFSfapprox(const MSSMNoFV_onshell* model)
{
   return gm2calc::amu2LFSfapprox(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

double gm2calc_amu2LFSfapprox_non_tan_beta_resummed(const MSSMNoFV_onshell* model)
{
   return gm2calc::amu2LFSfapprox_non_tan_beta_resummed(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

/* === photonic 2-loop corrections === */

double gm2calc_amuChipmPhotonic(const MSSMNoFV_onshell* model)
{
   return gm2calc::amuChipmPhotonic(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

double gm2calc_amuChi0Photonic(const MSSMNoFV_onshell* model)
{
   return gm2calc::amuChi0Photonic(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

/* === SUSY 2L(a) diagrams === */

double gm2calc_amu2LaSferm(const MSSMNoFV_onshell* model)
{
   return gm2calc::amu2LaSferm(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

double gm2calc_amu2LaCha(const MSSMNoFV_onshell* model)
{
   return gm2calc::amu2LaCha(
      *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model));
}

} /* extern "C" */
