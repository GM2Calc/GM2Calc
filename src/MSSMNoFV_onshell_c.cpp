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

#include "MSSMNoFV_onshell.h"
#include "MSSMNoFV_onshell.hpp"
#include "gm2_error.hpp"
#include <iostream>

extern "C"
{

const char* gm2calc_error_str(enum EError error)
{
   const char* error_str = "Unknown error";

   switch (error) {
   case NoError:
      error_str = "no error";
      break;
   case InvalidInput:
      error_str = "Input parameter set to invalid value";
      break;
   case PhysicalProblem:
      error_str = "Physical problem has occurred during calculation";
      break;
   default:
      break;
   }

   return error_str;
}

MSSMNoFV_onshell* gm2calc_mssmnofv_new()
{
   return reinterpret_cast<MSSMNoFV_onshell*>(new gm2calc::MSSMNoFV_onshell());
}

void gm2calc_mssmnofv_free(MSSMNoFV_onshell* model)
{
   delete reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model);
}

void gm2calc_mssmnofv_set_alpha_MZ(MSSMNoFV_onshell* model, double alpha_MZ)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_alpha_MZ(alpha_MZ);
}

void gm2calc_mssmnofv_set_alpha_thompson(MSSMNoFV_onshell* model, double alpha_0)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_alpha_thompson(alpha_0);
}

void gm2calc_mssmnofv_set_Ae(MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Ae(i,k,a);
}

void gm2calc_mssmnofv_set_Au(MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Au(i,k,a);
}

void gm2calc_mssmnofv_set_Ad(MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Ad(i,k,a);
}

void gm2calc_mssmnofv_set_g3(MSSMNoFV_onshell* model, double g3)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_g3(g3);
}

void gm2calc_mssmnofv_set_MassB(MSSMNoFV_onshell* model, double mass_b)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassB(mass_b);
}

void gm2calc_mssmnofv_set_MassWB(MSSMNoFV_onshell* model, double mass_wb)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassWB(mass_wb);
}

void gm2calc_mssmnofv_set_MassG(MSSMNoFV_onshell* model, double mass_g)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassG(mass_g);
}

void gm2calc_mssmnofv_set_mq2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double mq2)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_mq2(i,k,mq2);
}

void gm2calc_mssmnofv_set_mu2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double mu2)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_mu2(i,k,mu2);
}

void gm2calc_mssmnofv_set_md2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double md2)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_md2(i,k,md2);
}

void gm2calc_mssmnofv_set_ml2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double ml2)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_ml2(i,k,ml2);
}

void gm2calc_mssmnofv_set_me2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double me2)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_me2(i,k,me2);
}

void gm2calc_mssmnofv_set_Mu(MSSMNoFV_onshell* model, double mu)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Mu(mu);
}

void gm2calc_mssmnofv_set_TB(MSSMNoFV_onshell* model, double tan_beta)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_TB(tan_beta);
}

void gm2calc_mssmnofv_set_scale(MSSMNoFV_onshell* model, double scale)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_scale(scale);
}

void gm2calc_mssmnofv_set_MAh_pole(MSSMNoFV_onshell* model, double MA0)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MA0(MA0);
}

void gm2calc_mssmnofv_set_MZ_pole(MSSMNoFV_onshell* model, double MZ)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MVZ = MZ;
}

void gm2calc_mssmnofv_set_MW_pole(MSSMNoFV_onshell* model, double MW)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MVWm = MW;
}

void gm2calc_mssmnofv_set_MFt_pole(MSSMNoFV_onshell* model, double MFt)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MFt = MFt;
}

void gm2calc_mssmnofv_set_MFb_running(MSSMNoFV_onshell* model, double MFb)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MFb = MFb;
}

void gm2calc_mssmnofv_set_MFtau_pole(MSSMNoFV_onshell* model, double MFtau)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MFtau = MFtau;
}

void gm2calc_mssmnofv_set_MFm_pole(MSSMNoFV_onshell* model, double MFm)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MFm = MFm;
}

void gm2calc_mssmnofv_set_MSm_pole(MSSMNoFV_onshell* model, unsigned i, double MSm)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MSm(i) = MSm;
}

void gm2calc_mssmnofv_set_MSvmL_pole(MSSMNoFV_onshell* model, double MSvmL)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MSvmL = MSvmL;
}

void gm2calc_mssmnofv_set_MCha_pole(MSSMNoFV_onshell* model, unsigned i, double MCha)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MCha(i) = MCha;
}

void gm2calc_mssmnofv_set_MChi_pole(MSSMNoFV_onshell* model, unsigned i, double MChi)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->get_physical().MChi(i) = MChi;
}

void gm2calc_mssmnofv_set_verbose_output(MSSMNoFV_onshell* model, int verbose_output)
{
   reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_verbose_output(verbose_output);
}

int gm2calc_mssmnofv_convert_to_onshell(MSSMNoFV_onshell* model)
{
   int error = NoError;

   try {
      reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->convert_to_onshell();
   } catch (const gm2calc::EInvalidInput&) {
      error = InvalidInput;
   } catch (const gm2calc::EPhysicalProblem&) {
      error = InvalidInput;
   } catch (...) {
      error = UnknownError;
   }

   return error;
}

int gm2calc_mssmnofv_convert_to_onshell_params(
   MSSMNoFV_onshell* model, double precision, unsigned max_iterations)
{
   int error = NoError;

   try {
      reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->convert_to_onshell(precision, max_iterations);
   } catch (const gm2calc::EInvalidInput&) {
      error = InvalidInput;
   } catch (const gm2calc::EPhysicalProblem&) {
      error = InvalidInput;
   } catch (...) {
      error = UnknownError;
   }

   return error;
}

int gm2calc_mssmnofv_calculate_masses(MSSMNoFV_onshell* model)
{
   int error = NoError;

   try {
      reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->calculate_masses();
   } catch (const gm2calc::EInvalidInput&) {
      error = InvalidInput;
   } catch (const gm2calc::EPhysicalProblem&) {
      error = InvalidInput;
   } catch (...) {
      error = UnknownError;
   }

   return error;
}

void print_mssmnofv(const MSSMNoFV_onshell* model)
{
   std::cout << *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model);
}

} // extern "C"
