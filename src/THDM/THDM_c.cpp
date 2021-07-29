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

#include <complex>

namespace gm2calc {
namespace {

gm2calc::SM convert_to_SM(::SM* sm)
{
   gm2calc::SM s;

   if (sm != nullptr) {
      s.set_alpha_em_0(sm->alpha_em_0);
      s.set_alpha_em_mz(sm->alpha_em_mz);
      s.set_alpha_s_mz(sm->alpha_s_mz);
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
         s.set_mv(i, sm->mv[i]);
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

gm2calc::thdm::Gauge_basis convert_to_basis(THDM_gauge_basis* basis)
{
   gm2calc::thdm::Gauge_basis b;

   if (basis != nullptr) {
      switch (basis->yukawa_type) {
      case THDM_type_1:
         b.yukawa_type = thdm::Yukawa_type::type_1;
         break;
      case THDM_type_2:
         b.yukawa_type = thdm::Yukawa_type::type_2;
         break;
      case THDM_type_X:
         b.yukawa_type = thdm::Yukawa_type::type_X;
         break;
      case THDM_type_Y:
         b.yukawa_type = thdm::Yukawa_type::type_Y;
         break;
      case THDM_aligned:
         b.yukawa_type = thdm::Yukawa_type::aligned;
         break;
      case THDM_general:
         b.yukawa_type = thdm::Yukawa_type::general;
         break;
      }
      for (int i = 0; i < 7; i++) {
         b.lambda(i) = basis->lambda[i];
      }
      b.tan_beta = basis->tan_beta;
      b.m122 = basis->m122;
      b.zeta_u = basis->zeta_u;
      b.zeta_d = basis->zeta_d;
      b.zeta_l = basis->zeta_l;
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            b.Xu(i, k) = std::complex<double>(basis->Xu_real[i][k], basis->Xu_imag[i][k]);
            b.Xd(i, k) = std::complex<double>(basis->Xd_real[i][k], basis->Xd_imag[i][k]);
            b.Xl(i, k) = std::complex<double>(basis->Xl_real[i][k], basis->Xl_imag[i][k]);
         }
      }
   }

   return b;
}

gm2calc::thdm::Mass_basis convert_to_basis(THDM_mass_basis* basis)
{
   gm2calc::thdm::Mass_basis b;

   if (basis != nullptr) {
      switch (basis->yukawa_type) {
      case THDM_type_1:
         b.yukawa_type = thdm::Yukawa_type::type_1;
         break;
      case THDM_type_2:
         b.yukawa_type = thdm::Yukawa_type::type_2;
         break;
      case THDM_type_X:
         b.yukawa_type = thdm::Yukawa_type::type_X;
         break;
      case THDM_type_Y:
         b.yukawa_type = thdm::Yukawa_type::type_Y;
         break;
      case THDM_aligned:
         b.yukawa_type = thdm::Yukawa_type::aligned;
         break;
      case THDM_general:
         b.yukawa_type = thdm::Yukawa_type::general;
         break;
      }
      b.mh = basis->mh;
      b.mH = basis->mH;
      b.mA = basis->mA;
      b.mHp = basis->mHp;
      b.sin_beta_minus_alpha = basis->sin_beta_minus_alpha;
      b.lambda_6 = basis->lambda_6;
      b.lambda_7 = basis->lambda_7;
      b.tan_beta = basis->tan_beta;
      b.m122 = basis->m122;
      b.zeta_u = basis->zeta_u;
      b.zeta_d = basis->zeta_d;
      b.zeta_l = basis->zeta_l;
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            b.Xu(i, k) = std::complex<double>(basis->Xu_real[i][k], basis->Xu_imag[i][k]);
            b.Xd(i, k) = std::complex<double>(basis->Xd_real[i][k], basis->Xd_imag[i][k]);
            b.Xl(i, k) = std::complex<double>(basis->Xl_real[i][k], basis->Xl_imag[i][k]);
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
gm2calc_error gm2calc_thdm_new_with_gauge_basis(THDM** model, THDM_gauge_basis* basis, ::SM* sm)
{
   if (model == nullptr) {
      return gm2calc_InvalidInput;
   }

   gm2calc_error error;

   try {
      const auto b(gm2calc::convert_to_basis(basis));
      const auto s(gm2calc::convert_to_SM(sm));
      *model = reinterpret_cast<THDM*>(new gm2calc::THDM(b, s));
      error = gm2calc_NoError;
   } catch (const gm2calc::EInvalidInput&) {
      *model = nullptr;
      error = gm2calc_InvalidInput;
   } catch (const gm2calc::EPhysicalProblem&) {
      *model = nullptr;
      error = gm2calc_PhysicalProblem;
   } catch (...) {
      *model = nullptr;
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
gm2calc_error gm2calc_thdm_new_with_mass_basis(THDM** model, THDM_mass_basis* basis, ::SM* sm)
{
   if (model == nullptr) {
      return gm2calc_InvalidInput;
   }

   gm2calc_error error;

   try {
      const auto b(gm2calc::convert_to_basis(basis));
      const auto s(gm2calc::convert_to_SM(sm));
      *model = reinterpret_cast<THDM*>(new gm2calc::THDM(b, s));
      error = gm2calc_NoError;
   } catch (const gm2calc::EInvalidInput&) {
      *model = nullptr;
      error = gm2calc_InvalidInput;
   } catch (const gm2calc::EPhysicalProblem&) {
      *model = nullptr;
      error = gm2calc_PhysicalProblem;
   } catch (...) {
      *model = nullptr;
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
