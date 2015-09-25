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

#include "gm2_1loop.hpp"
#include "gm2_2loop.hpp"
#include "gm2_error.hpp"
#include "MSSMNoFV_onshell.hpp"

#include <cstdio>
#include <iostream>

double calculate_amu(gm2calc::MSSMNoFV_onshell& model)
{
   return gm2calc::calculate_amu_1loop(model)
      + gm2calc::calculate_amu_2loop(model);
}

gm2calc::MSSMNoFV_onshell setup()
{
   gm2calc::MSSMNoFV_onshell model;

   const double g3 = 1.06459;
   const double MW = 80.385;
   const double MZ = 91.1876;
   const double ME = 0.00051;
   const double MM = 0.105658;
   const double ML = 1.777;
   const double MT = 173.5;
   const double MB = 4.18;
   Eigen::Matrix<double,3,3> Ae, Au, Ad;
   Eigen::Matrix<double,3,3> mq2, ml2, md2, mu2, me2;
   const double Mu = 350;
   const double M1 = 150;
   const double M2 = 300;
   const double M3 = 1;
   const double scale = 454.7;
   const double MA0 = 1500;
   const double TB = 10;

   Au.setZero();
   Ad.setZero();
   Ae.setZero();

   mq2.diagonal().setConstant(500 * 500);
   ml2 = md2 = mu2 = me2 = mq2;

   // set parameters
   model.get_physical().MVWm = MW;
   model.get_physical().MVZ = MZ;
   model.get_physical().MFe = ME;
   model.get_physical().MFm = MM;
   model.get_physical().MFtau = ML;
   model.get_physical().MFt = MT;
   model.get_physical().MFb = MB;

   model.set_g3(g3);
   model.set_TB(TB);
   model.set_Ae(Ae);
   model.set_Ad(Ad);
   model.set_Au(Au);
   model.set_Mu(Mu);
   model.set_mq2(mq2);
   model.set_ml2(ml2);
   model.set_md2(md2);
   model.set_mu2(mu2);
   model.set_me2(me2);
   model.set_MassB(M1);
   model.set_MassWB(M2);
   model.set_MassG(M3);
   model.set_scale(scale);
   model.set_MA0(MA0);

   return model;
}

int main()
{
   const double tanb_start = 2.;
   const double tanb_stop = 100.;
   const unsigned nsteps = 100;

   printf("# %14s %16s %16s\n", "tan(beta)", "amu", "error");

   for (unsigned n = 0; n < nsteps; n++) {
      double amu;
      const double tanb = tanb_start + (tanb_stop - tanb_start) * n / nsteps;
      std::string error;

      gm2calc::MSSMNoFV_onshell model(setup());
      model.do_force_output(false); // throw exception in case of problem
      model.set_TB(tanb);

      try {
         model.calculate_masses();
         amu = calculate_amu(model);
      } catch (const gm2calc::Error& e) {
         error = "# " + e.what();
         amu = std::numeric_limits<double>::signaling_NaN();
      }

      printf("%16.8e %16.8e %s\n", tanb, amu, error.c_str());
   }

   return 0;
}
