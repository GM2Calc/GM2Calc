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
#include <iostream>

extern "C"
{

MSSMNoFV_onshell* gm2calc_mssmnofv_new()
{
   return reinterpret_cast<MSSMNoFV_onshell*>(new gm2calc::MSSMNoFV_onshell());
}

void gm2calc_mssmnofv_free(MSSMNoFV_onshell* model)
{
   delete reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model);
}

void gm2calc_mssmnofv_set_TB(MSSMNoFV_onshell* model, double tan_beta)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_TB(tan_beta);
}

void gm2calc_mssmnofv_set_Ae(MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Ae(i,k,a);
}

void gm2calc_mssmnofv_set_Au(MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Au(i,k,a);
}

void gm2calc_mssmnofv_set_Ad(MSSMNoFV_onshell* model, unsigned i, unsigned k, double a)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Ad(i,k,a);
}

void gm2calc_mssmnofv_set_Mu(MSSMNoFV_onshell* model, double mu)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_Mu(mu);
}

void gm2calc_mssmnofv_set_MassB(MSSMNoFV_onshell* model, double mass_b)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassB(mass_b);
}

void gm2calc_mssmnofv_set_MassWB(MSSMNoFV_onshell* model, double mass_wb)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassWB(mass_wb);
}

void gm2calc_mssmnofv_set_MassG(MSSMNoFV_onshell* model, double mass_g)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MassG(mass_g);
}

void gm2calc_mssmnofv_set_MA0(MSSMNoFV_onshell* model, double MA0)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_MA0(MA0);
}

void gm2calc_mssmnofv_set_mq2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double mq2)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_mq2(i,k,mq2);
}

void gm2calc_mssmnofv_set_mu2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double mu2)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_mu2(i,k,mu2);
}

void gm2calc_mssmnofv_set_md2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double md2)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_md2(i,k,md2);
}

void gm2calc_mssmnofv_set_ml2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double ml2)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_ml2(i,k,ml2);
}

void gm2calc_mssmnofv_set_me2(MSSMNoFV_onshell* model, unsigned i, unsigned k, double me2)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_me2(i,k,me2);
}

void gm2calc_mssmnofv_set_scale(MSSMNoFV_onshell* model, double scale)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->set_scale(scale);
}

void gm2calc_mssmnofv_calculate_masses(MSSMNoFV_onshell* model)
{
   return reinterpret_cast<gm2calc::MSSMNoFV_onshell*>(model)->calculate_masses();
}

void print_mssmnofv(const MSSMNoFV_onshell* model)
{
   std::cout << *reinterpret_cast<const gm2calc::MSSMNoFV_onshell*>(model);
}

} // extern "C"
