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

#include "gm2calc/THDM.h"
#include "gm2calc/THDM.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2calc/SM.h"

#include "gm2_log.hpp"

#include <complex>
#include <string>
#include <iostream>

namespace gm2calc {
namespace {

gm2calc::SM convert_to_SM(::SM* sm)
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

gm2calc::THDM::General_basis convert_to_basis(THDM_general_basis* basis)
{
   gm2calc::THDM::General_basis b;

   if (basis != 0) {
      switch (basis->yukawa_scheme) {
      case THDM_type_1:
         b.yukawa_scheme = THDM::Yukawa_scheme::type_1;
         break;
      case THDM_type_2:
         b.yukawa_scheme = THDM::Yukawa_scheme::type_2;
         break;
      case THDM_type_X:
         b.yukawa_scheme = THDM::Yukawa_scheme::type_X;
         break;
      case THDM_type_Y:
         b.yukawa_scheme = THDM::Yukawa_scheme::type_Y;
         break;
      case THDM_general:
         b.yukawa_scheme = THDM::Yukawa_scheme::general;
         break;
      }
      b.lambda1 = basis->lambda1;
      b.lambda2 = basis->lambda2;
      b.lambda3 = basis->lambda3;
      b.lambda4 = basis->lambda4;
      b.lambda5 = basis->lambda5;
      b.lambda6 = basis->lambda6;
      b.lambda7 = basis->lambda7;
      b.tan_beta = basis->tan_beta;
      b.m122 = basis->m122;
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            b.Xu(i, k) = basis->Xu[i][k];
            b.Xd(i, k) = basis->Xd[i][k];
            b.Xl(i, k) = basis->Xl[i][k];
         }
      }
   }

   return b;
}

gm2calc::THDM::Physical_basis convert_to_basis(THDM_physical_basis* basis)
{
   gm2calc::THDM::Physical_basis b;

   if (basis != 0) {
      switch (basis->yukawa_scheme) {
      case THDM_type_1:
         b.yukawa_scheme = THDM::Yukawa_scheme::type_1;
         break;
      case THDM_type_2:
         b.yukawa_scheme = THDM::Yukawa_scheme::type_2;
         break;
      case THDM_type_X:
         b.yukawa_scheme = THDM::Yukawa_scheme::type_X;
         break;
      case THDM_type_Y:
         b.yukawa_scheme = THDM::Yukawa_scheme::type_Y;
         break;
      case THDM_general:
         b.yukawa_scheme = THDM::Yukawa_scheme::general;
         break;
      }
      b.mh = basis->mh;
      b.mH = basis->mH;
      b.mA = basis->mA;
      b.mHp = basis->mHp;
      b.sin_beta_minus_alpha = basis->sin_beta_minus_alpha;
      b.lambda6 = basis->lambda6;
      b.lambda7 = basis->lambda7;
      b.tan_beta = basis->tan_beta;
      b.m122 = basis->m122;
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            b.Xu(i, k) = basis->Xu[i][k];
            b.Xd(i, k) = basis->Xd[i][k];
            b.Xl(i, k) = basis->Xl[i][k];
         }
      }
   }

   return b;
}

} // anonymous namespace
} // namespace gm2calc

/**
 * @file THDM_c.cpp
 * @brief contains definitions of C interface functions for the model
 *
 * This file contains the definitions for the C interface functions
 * used to set and retrieve the model parameters and masses.
 */

extern "C"
{

/**
 * @brief Allocate a new general THDM model with general basis input.
 *
 * This function allocates a new general THDM model and sets the
 * pointer to the created object.  To prevent a resource leak, the
 * model should be destroyed using gm2calc_thdm_free() .
 *
 * If an error occurs, the model pointer will be set to 0 and a
 * corresponding error code will be returned.
 *
 * @return error code
 */
gm2calc_error gm2calc_thdm_new_with_general_basis(THDM** model, THDM_general_basis* basis, ::SM* sm)
{
   if (model == 0) {
      return gm2calc_InvalidInput;
   }

   gm2calc_error error;

   try {
      const auto b(gm2calc::convert_to_basis(basis));
      const auto s(gm2calc::convert_to_SM(sm));
      *model = reinterpret_cast<THDM*>(new gm2calc::THDM(b, s));
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
 * model should be destroyed using gm2calc_thdm_free() .
 *
 * If an error occurs, the model pointer will be set to 0 and a
 * corresponding error code will be returned.
 *
 * @return error code
 */
gm2calc_error gm2calc_thdm_new_with_physical_basis(THDM** model, THDM_physical_basis* basis, ::SM* sm)
{
   if (model == 0) {
      return gm2calc_InvalidInput;
   }

   gm2calc_error error;

   try {
      const auto b(gm2calc::convert_to_basis(basis));
      const auto s(gm2calc::convert_to_SM(sm));
      *model = reinterpret_cast<THDM*>(new gm2calc::THDM(b, s));
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
 * created using gm2calc_thdm_new() .
 *
 * @param model pointer to model object
 */
void gm2calc_thdm_free(THDM* model)
{
   delete reinterpret_cast<gm2calc::THDM*>(model);
}

} // extern "C"
