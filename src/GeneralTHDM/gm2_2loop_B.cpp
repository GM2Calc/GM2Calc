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

#include "GeneralTHDM/gm2_2loop_helpers.hpp"

/**
 * \file gm2_2loop_B.cpp
 *
 * Contains functions necessary to calculate the bosonic THDM
 * contributions for g-2 at the 2-loop level.
 */

namespace gm2calc {

namespace general_thdm {

/**
 * Calculates 2-loop bosonic pure electroweak contributions.
 *
 * Eq (49), arxiv:1607:06292
 */
double amu2L_B_EWadd(double eta, double zetal)
{
   return 2.3e-11 * eta * zetal;
}

/**
 * Calculates 2-loop bosonic non-Yukawa contributions.
 *
 * Eq (71), arxiv:1607:06292
 */
double amu2L_B_nonYuk()
{
   // @todo(alex) implementation missing
   return 0.0;
}

/**
 * Calculates 2-loop bosonic Yukawa contributions.
 *
 * Eq (52), arxiv:1607:06292
 */
double amu2L_B_Yuk()
{
   // @todo(alex) implementation missing

   const auto tb = 1.0;
   const auto zetal = 0.0;
   const auto lambda5 = 0.0;
   const auto eta = 0.0;

   const auto a000 = 0.0;
   const auto a0z0 = 0.0;
   const auto a500 = 0.0;
   const auto a5z0 = 0.0;
   const auto a001 = 0.0;
   const auto a0z1 = 0.0;
   const auto a501 = 0.0;
   const auto a5z1 = 0.0;

   const double res =
      + a000
      + a0z0*(tb - 1.0/tb)*zetal
      + a500*lambda5
      + a5z0*(tb - 1.0/tb)*lambda5*zetal
      + (
         + a001*(tb - 1.0/tb)
         + a0z1*zetal
         + a501*(tb - 1.0/tb)*lambda5
         + a5z1*lambda5*zetal
         )*eta;

   return res;
}

} // namespace general_thdm

} // namespace gm2calc
