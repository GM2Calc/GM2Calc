#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2calc/GeneralTHDM.hpp"

#include <iostream>

int main()
{
   gm2calc::GeneralTHDM::Physical_basis basis;
   basis.yukawa_scheme = gm2calc::GeneralTHDM::Yukawa_scheme::type_2;
   basis.mh = 125;
   basis.mH = 400;
   basis.mA = 420;
   basis.mHp = 440;
   basis.sin_beta_minus_alpha = 0.999;
   basis.lambda6 = 0;
   basis.lambda7 = 0;
   basis.tan_beta = 3;
   basis.m122 = 40000;

   gm2calc::SM sm;
   sm.set_alpha_em_mz(1.0/128.94579);
   sm.set_mu(2, 173.34);
   sm.set_mu(1, 1.28);
   sm.set_md(2, 4.18);
   sm.set_ml(2, 1.77684);

   try {
      gm2calc::GeneralTHDM model(basis, sm);

      const double amu = gm2calc::calculate_amu_1loop(model)
                       + gm2calc::calculate_amu_2loop(model);

      const double delta_amu =
         gm2calc::calculate_uncertainty_amu_2loop(model);

      std::cout << "amu = " << amu << " +- " << delta_amu << '\n';
   } catch (const gm2calc::Error& e) {
      std::cout << e.what() << '\n';
   }

   return 0;
}