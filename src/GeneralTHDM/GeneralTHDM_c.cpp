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
#include "gm2calc/SM.h"

#include "gm2_log.hpp"

#include <complex>
#include <string>
#include <iostream>

namespace gm2calc {
namespace {

gm2calc::SM convert_to_SM(gm2calc_SM* sm)
{
   gm2calc::SM s;

   if (sm != 0) {
      s.set_alpha_em_0(sm->alpha_em_0);
      s.set_alpha_em_mz(sm->alpha_em_mz);
      s.set_mh(sm->mh);
      s.set_mw(sm->mw);
      s.set_mz(sm->mz);
      for (int i = 0; i < 3; i++) {
         s.set_mu(i, sm->mu[i]);
      }
      for (int i = 0; i < 3; i++) {
         s.set_md(i, sm->md[i]);
      }
      for (int i = 0; i < 3; i++) {
         s.set_ml(i, sm->ml[i]);
      }
      Eigen::Matrix<std::complex<double>,3,3> ckm;
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            ckm(i, k) = std::complex<double>(sm->ckm_real[i][k], sm->ckm_imag[i][k]);
         }
      }
      s.set_ckm(ckm);
   }

   return s;
}

} // anonymous namespace
} // namespace gm2calc

/**
 * @file GeneralTHDM_c.cpp
 * @brief contains definitions of C interface functions for the model
 *
 * This file contains the definitions for the C interface functions
 * used to set and retrieve the model parameters and masses.
 */

extern "C"
{

/* ********** general basis ********** */

/**
 * @brief Allocate a new general basis for the general THDM.
 * @return pointer to general basis
 */
GeneralTHDM_general_basis* gm2calc_generalthdm_general_basis_new()
{
   return reinterpret_cast<GeneralTHDM_general_basis*>(new gm2calc::GeneralTHDM::General_basis());
}

/**
 * @brief Deletes a general basis for the general THDM.
 * @param pointer to general basis
 */
void gm2calc_generalthdm_general_basis_free(GeneralTHDM_general_basis* basis)
{
   delete reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis);
}

void gm2calc_generalthdm_general_basis_set_yukawa_scheme(GeneralTHDM_general_basis* basis, GeneralTHDM_yukawa_scheme yukawa_scheme)
{
   auto b = reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis);

   switch (yukawa_scheme) {
   case GeneralTHDM_type_1:
      b->yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::type_1;
      break;
   case GeneralTHDM_type_2:
      b->yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::type_2;
      break;
   case GeneralTHDM_type_X:
      b->yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::type_X;
      break;
   case GeneralTHDM_type_Y:
      b->yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::type_Y;
      break;
   case GeneralTHDM_general:
      b->yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::general;
      break;
   }
}

void gm2calc_generalthdm_general_basis_set_lambda_1(GeneralTHDM_general_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->lambda1 = val;
}

void gm2calc_generalthdm_general_basis_set_lambda_2(GeneralTHDM_general_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->lambda2 = val;
}

void gm2calc_generalthdm_general_basis_set_lambda_3(GeneralTHDM_general_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->lambda3 = val;
}

void gm2calc_generalthdm_general_basis_set_lambda_4(GeneralTHDM_general_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->lambda4 = val;
}

void gm2calc_generalthdm_general_basis_set_lambda_5(GeneralTHDM_general_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->lambda5 = val;
}

void gm2calc_generalthdm_general_basis_set_lambda_6(GeneralTHDM_general_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->lambda6 = val;
}

void gm2calc_generalthdm_general_basis_set_lambda_7(GeneralTHDM_general_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->lambda7 = val;
}

void gm2calc_generalthdm_general_basis_set_tan_beta(GeneralTHDM_general_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->tan_beta = val;
}

void gm2calc_generalthdm_general_basis_set_m122(GeneralTHDM_general_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->m122 = val;
}

void gm2calc_generalthdm_general_basis_set_Xu(GeneralTHDM_general_basis* basis, unsigned i, unsigned j, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->Xu(i, j) = val;
}

void gm2calc_generalthdm_general_basis_set_Xd(GeneralTHDM_general_basis* basis, unsigned i, unsigned j, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->Xd(i, j) = val;
}

void gm2calc_generalthdm_general_basis_set_Xl(GeneralTHDM_general_basis* basis, unsigned i, unsigned j, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::General_basis*>(basis)->Xl(i, j) = val;
}

/* ********** physical basis ********** */

/**
 * @brief Allocate a new physical basis for the general THDM.
 * @return pointer to physical basis
 */
GeneralTHDM_physical_basis* gm2calc_generalthdm_physical_basis_new()
{
   return reinterpret_cast<GeneralTHDM_physical_basis*>(new gm2calc::GeneralTHDM::Physical_basis());
}

/**
 * @brief Deletes a physical basis for the physical THDM.
 * @param pointer to physical basis
 */
void gm2calc_generalthdm_physical_basis_free(GeneralTHDM_physical_basis* basis)
{
   delete reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis);
}

void gm2calc_generalthdm_physical_basis_set_yukawa_scheme(GeneralTHDM_physical_basis* basis, GeneralTHDM_yukawa_scheme yukawa_scheme)
{
   auto b = reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis);

   switch (yukawa_scheme) {
   case GeneralTHDM_type_1:
      b->yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::type_1;
      break;
   case GeneralTHDM_type_2:
      b->yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::type_2;
      break;
   case GeneralTHDM_type_X:
      b->yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::type_X;
      break;
   case GeneralTHDM_type_Y:
      b->yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::type_Y;
      break;
   case GeneralTHDM_general:
      b->yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::general;
      break;
   }
}

void gm2calc_generalthdm_physical_basis_set_mh(GeneralTHDM_physical_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->mh = val;
}

void gm2calc_generalthdm_physical_basis_set_mH(GeneralTHDM_physical_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->mH = val;
}

void gm2calc_generalthdm_physical_basis_set_mA(GeneralTHDM_physical_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->mA = val;
}

void gm2calc_generalthdm_physical_basis_set_mHp(GeneralTHDM_physical_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->mHp = val;
}

void gm2calc_generalthdm_physical_basis_set_sin_beta_minus_alpha(GeneralTHDM_physical_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->sin_beta_minus_alpha = val;
}

void gm2calc_generalthdm_physical_basis_set_lambda_6(GeneralTHDM_physical_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->lambda6 = val;
}

void gm2calc_generalthdm_physical_basis_set_lambda_7(GeneralTHDM_physical_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->lambda7 = val;
}

void gm2calc_generalthdm_physical_basis_set_tan_beta(GeneralTHDM_physical_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->tan_beta = val;
}

void gm2calc_generalthdm_physical_basis_set_m122(GeneralTHDM_physical_basis* basis, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->m122 = val;
}

void gm2calc_generalthdm_physical_basis_set_Xu(GeneralTHDM_physical_basis* basis, unsigned i, unsigned j, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->Xu(i, j) = val;
}

void gm2calc_generalthdm_physical_basis_set_Xd(GeneralTHDM_physical_basis* basis, unsigned i, unsigned j, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->Xd(i, j) = val;
}

void gm2calc_generalthdm_physical_basis_set_Xl(GeneralTHDM_physical_basis* basis, unsigned i, unsigned j, double val)
{
   reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis)->Xl(i, j) = val;
}

/* ********** THDM model ********** */

/**
 * @brief Allocate a new general THDM model with general basis input.
 *
 * This function allocates a new general THDM model and sets the
 * pointer to the created object.  To prevent a resource leak, the
 * model should be destroyed using gm2calc_generalthdm_free() .
 *
 * If an error occurs, the model pointer will be set to 0 and a
 * corresponding error code will be returned.
 *
 * @return error code
 */
gm2calc_error gm2calc_generalthdm_new_with_general_basis(GeneralTHDM** model, GeneralTHDM_general_basis* basis, gm2calc_SM* sm)
{
   if (model == 0 || basis == 0 || sm == 0) {
      return gm2calc_InvalidInput;
   }

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
 * If an error occurs, the model pointer will be set to 0 and a
 * corresponding error code will be returned.
 *
 * @return error code
 */
gm2calc_error gm2calc_generalthdm_new_with_physical_basis(GeneralTHDM** model, GeneralTHDM_physical_basis* basis, gm2calc_SM* sm)
{
   if (model == 0 || basis == 0 || sm == 0) {
      return gm2calc_InvalidInput;
   }

   gm2calc_error error;

   try {
      const auto b = reinterpret_cast<gm2calc::GeneralTHDM::Physical_basis*>(basis);
      const gm2calc::SM s(gm2calc::convert_to_SM(sm));
      *model = reinterpret_cast<GeneralTHDM*>(new gm2calc::GeneralTHDM(*b, s));
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
