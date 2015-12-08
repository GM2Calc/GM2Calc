#include "gm2_1loop.hpp"
#include "gm2_2loop.hpp"
#include "MSSMNoFV_onshell.hpp"
#include <iostream>

gm2calc::MSSMNoFV_onshell setup() {
   gm2calc::MSSMNoFV_onshell model;

   const Eigen::Matrix<double,3,3> UnitMatrix
      = Eigen::Matrix<double,3,3>::Identity();

   // fill DR-bar parameters
   model.set_TB(10);                      // 1L
   model.set_Ae(1,1,0);                   // 1L

   // fill on-shell parameters
   model.set_Mu(350);                     // 1L
   model.set_MassB(150);                  // 1L
   model.set_MassWB(300);                 // 1L
   model.set_MassG(1000);                 // 2L
   model.set_mq2(500 * 500 * UnitMatrix); // 2L
   model.set_ml2(500 * 500 * UnitMatrix); // 1L(smuon)/2L
   model.set_md2(500 * 500 * UnitMatrix); // 2L
   model.set_mu2(500 * 500 * UnitMatrix); // 2L
   model.set_me2(500 * 500 * UnitMatrix); // 1L(smuon)/2L
   model.set_Au(2,2,0);                   // 2L
   model.set_Ad(2,2,0);                   // 2L
   model.set_Ae(2,2,0);                   // 2L
   model.set_MA0(1500);                   // 2L
   model.set_scale(454.7);                // 2L

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
