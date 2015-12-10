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

extern "C"
{

C_MSSMNoFV_onshell* c_mssm_onshell_new()
{
   return reinterpret_cast<C_MSSMNoFV_onshell*>(new gm2calc::MSSMNoFV_onshell());
}

void c_mssm_onshell_free(C_MSSMNoFV_onshell* model)
{
   delete reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model);
}

void c_mssm_onshell_set_TB(C_MSSMNoFV_onshell* model, double tan_beta)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_TB(tan_beta);
}

void c_mssm_onshell_set_Ae(C_MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Ae(i,k,a);
}

void c_mssm_onshell_set_Au(C_MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Au(i,k,a);
}

void c_mssm_onshell_set_Ad(C_MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Ad(i,k,a);
}

void c_mssm_onshell_set_Mu(C_MSSMNoFV_onshell* model, double mu)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Mu(mu);
}

void c_mssm_onshell_set_MassB(C_MSSMNoFV_onshell* model, double mass_b)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassB(mass_b);
}

void c_mssm_onshell_set_MassWB(C_MSSMNoFV_onshell* model, double mass_wb)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassWB(mass_wb);
}

void c_mssm_onshell_set_MassG(C_MSSMNoFV_onshell* model, double mass_g)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassG(mass_g);
}

void c_mssm_onshell_set_MA0(C_MSSMNoFV_onshell* model, double MA0)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MA0(MA0);
}

void c_mssm_onshell_set_mq2(C_MSSMNoFV_onshell* model, unsigned i, unsigned k, double mq2)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_mq2(i,k,mq2);
}

void c_mssm_onshell_set_mu2(C_MSSMNoFV_onshell* model, unsigned i, unsigned k, double mu2)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_mu2(i,k,mu2);
}

void c_mssm_onshell_set_md2(C_MSSMNoFV_onshell* model, unsigned i, unsigned k, double md2)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_md2(i,k,md2);
}

void c_mssm_onshell_set_ml2(C_MSSMNoFV_onshell* model, unsigned i, unsigned k, double ml2)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_ml2(i,k,ml2);
}

void c_mssm_onshell_set_me2(C_MSSMNoFV_onshell* model, unsigned i, unsigned k, double me2)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_me2(i,k,me2);
}

void c_mssm_onshell_set_scale(C_MSSMNoFV_onshell* model, double scale)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_scale(scale);
}

void c_mssm_onshell_calculate_masses(C_MSSMNoFV_onshell* model)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->calculate_masses();
}

} // extern "C"
