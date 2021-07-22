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

} // extern "C"
