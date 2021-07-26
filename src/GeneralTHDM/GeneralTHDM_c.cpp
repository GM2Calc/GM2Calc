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

#include "gm2calc/GeneralTHDM.h"
#include "gm2calc/GeneralTHDM.hpp"
#include "gm2calc/gm2_error.hpp"

#include "gm2_log.hpp"

#include <complex>
#include <string>
#include <iostream>

/**
 * @file GeneralTHDM_c.cpp
 * @brief contains definitions of C interface functions for the model
 *
 * This file contains the definitions for the C interface functions
 * used to set and retrieve the model parameters and masses.
 */

extern "C"
{

/**
 * @brief Allocate a new general basis for the general THDM.
 * @return pointer to general basis
 */
GeneralTHDM_general_basis* gm2calc_generalthdm_general_basis_new()
{
   return reinterpret_cast<GeneralTHDM_general_basis*>(new gm2calc::GeneralTHDM::General_basis());
}

/**
 * @brief Allocate a new physical basis for the general THDM.
 * @return pointer to physical basis
 */
GeneralTHDM_physical_basis* gm2calc_generalthdm_physical_basis_new()
{
   return reinterpret_cast<GeneralTHDM_physical_basis*>(new gm2calc::GeneralTHDM::Physical_basis());
}

/**
 * @brief Allocate a new general THDM model with general basis input.
 *
 * This function allocates a new general THDM model and sets the
 * pointer to the created object.  To prevent a resource leak, the
 * model should be destroyed using gm2calc_generalthdm_free() .
 *
 * @return error code
 */
gm2calc_error gm2calc_generalthdm_new_with_general_basis(GeneralTHDM** model, GeneralTHDM_general_basis* basis, gm2calc_SM* sm)
{
   gm2calc_error error;

   try {
      const auto b = reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis);
      const auto s = reinterpret_cast<gm2calc::SM*>(sm);
      *model = reinterpret_cast<GeneralTHDM*>(new gm2calc::GeneralTHDM(*b, *s));
      error = gm2calc_NoError;
   } catch (const gm2calc::EInvalidInput&) {
      *model = 0;
      error = gm2calc_InvalidInput;
   } catch (const gm2calc::EPhysicalProblem&) {
      *model = 0;
      error = gm2calc_PhysicalProblem;
   } catch (...) {
      *model = 0;
      error = gm2calc_UnknownError;
   }

   return error;
}

/**
 * @brief Allocate a new general THDM model with physical basis input.
 *
 * This function allocates a new general THDM model and sets the
 * pointer to the created object.  To prevent a resource leak, the
 * model should be destroyed using gm2calc_generalthdm_free() .
 *
 * @return error code
 */
gm2calc_error gm2calc_generalthdm_new_with_physical_basis(GeneralTHDM** model, GeneralTHDM_physical_basis* basis, gm2calc_SM* sm)
{
   gm2calc_error error;

   try {
      const auto b = reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis);
      const auto s = reinterpret_cast<gm2calc::SM*>(sm);
      *model = reinterpret_cast<GeneralTHDM*>(new gm2calc::GeneralTHDM(*b, *s));
      error = gm2calc_NoError;
   } catch (const gm2calc::EInvalidInput&) {
      *model = 0;
      error = gm2calc_InvalidInput;
   } catch (const gm2calc::EPhysicalProblem&) {
      *model = 0;
      error = gm2calc_PhysicalProblem;
   } catch (...) {
      *model = 0;
      error = gm2calc_UnknownError;
   }

   return error;
}

/**
 * @brief Deletes a general basis for the general THDM.
 * @param pointer to general basis
 */
void gm2calc_generalthdm_general_basis_free(GeneralTHDM_general_basis* basis)
{
   delete reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis);
}

/**
 * @brief Deletes a physical basis for the physical THDM.
 * @param pointer to physical basis
 */
void gm2calc_generalthdm_physical_basis_free(GeneralTHDM_physical_basis* basis)
{
   delete reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis);
}

/**
 * @brief Deletes a general THDM model.
 *
 * This function deletes a general THDM model object, which has been
 * created using gm2calc_generalthdm_new() .
 *
 * @param model pointer to model object
 */
void gm2calc_generalthdm_free(GeneralTHDM* model)
{
   delete reinterpret_cast<gm2calc::GeneralTHDM*>(model);
}

} // extern "C"
