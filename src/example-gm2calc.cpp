#include "gm2_1loop.hpp"
#include "gm2_2loop.hpp"
#include "MSSMNoFV_onshell.hpp"
#include <iostream>

gm2calc::MSSMNoFV_onshell setup() {
   gm2calc::MSSMNoFV_onshell model;

   const Eigen::Matrix<double,3,3> UnitMatrix
      = Eigen::Matrix<double,3,3>::Identity();

   // fill DR-bar parameters
   model.set_TB(10);
   model.set_Ae(1,1,0);

   // fill on-shell parameters
   model.set_Mu(350);
   model.set_MassB(150);
   model.set_MassWB(300);
   model.set_MassG(1000);
   model.set_mq2(500 * 500 * UnitMatrix);
   model.set_ml2(500 * 500 * UnitMatrix);
   model.set_md2(500 * 500 * UnitMatrix);
   model.set_mu2(500 * 500 * UnitMatrix);
   model.set_me2(500 * 500 * UnitMatrix);
   model.set_Au(2,2,0);
   model.set_Ad(2,2,0);
   model.set_Ae(2,2,0);
   model.set_MA0(1500);
   model.set_scale(454.7);

   // calculate mass spectrum
   model.calculate_masses();

   return model;
}

int main() {
   gm2calc::MSSMNoFV_onshell model(setup());

   const double amu =
      + gm2calc::calculate_amu_1loop(model)
      + gm2calc::calculate_amu_2loop(model);

   std::cout << "amu = " << amu << std::endl;

   return 0;
}
