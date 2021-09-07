#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"

#include "gm2calc/gm2_1loop.h"
#include "gm2calc/gm2_2loop.h"
#include "gm2calc/gm2_uncertainty.h"
#include "gm2calc/THDM.h"
#include "gm2calc/SM.h"

#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/THDM.hpp"

#include <utility>


void setup_SM(gm2calc::SM& cppsm, SM& csm)
{
   cppsm.set_alpha_em_0(1.0/137);
   cppsm.set_alpha_em_mz(1.0/130);
   cppsm.set_alpha_s_mz(0.12);
   cppsm.set_mh(125);
   cppsm.set_mw(80);
   cppsm.set_mz(91);
   cppsm.set_mu(2, 173);
   cppsm.set_mu(1, 1.3);
   cppsm.set_md(2, 4.2);
   cppsm.set_ml(2, 1.8);
   cppsm.set_ckm_from_wolfenstein(0.8, 0.7, 0.6, 0.5);

   csm.alpha_em_0 = cppsm.get_alpha_em_0();
   csm.alpha_em_mz = cppsm.get_alpha_em_mz();
   csm.alpha_s_mz = cppsm.get_alpha_s_mz();
   csm.mh = cppsm.get_mh();
   csm.mw = cppsm.get_mw();
   csm.mz = cppsm.get_mz();
   for (int i = 0; i < 3; i++) {
      csm.mu[i] = cppsm.get_mu(i);
      csm.md[i] = cppsm.get_md(i);
      csm.ml[i] = cppsm.get_ml(i);
   }
   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         csm.ckm_real[i][k] = std::real(cppsm.get_ckm(i, k));
         csm.ckm_imag[i][k] = std::imag(cppsm.get_ckm(i, k));
      }
   }
}


std::pair<gm2calc::THDM, THDM*> setup_mass_basis(gm2calc::thdm::Yukawa_type yukawa_type)
{
   gm2calc::thdm::Mass_basis basis;
   basis.yukawa_type = yukawa_type;
   basis.mh = 125;
   basis.mH = 400;
   basis.mA = 420;
   basis.mHp = 440;
   basis.sin_beta_minus_alpha = 0.999;
   basis.lambda_6 = 0;
   basis.lambda_7 = 0;
   basis.tan_beta = 3;
   basis.m122 = 40000;
   basis.zeta_u = 0.1;
   basis.zeta_d = 0.2;
   basis.zeta_l = 0.3;
   basis.Delta_u << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9;
   basis.Delta_d = 2.0*basis.Delta_u;
   basis.Delta_l = 3.0*basis.Delta_u;
   basis.Pi_u << 4.0*basis.Delta_u;
   basis.Pi_d = 2.0*basis.Pi_u;
   basis.Pi_l = 3.0*basis.Pi_u;

   THDM_mass_basis cbasis;
   cbasis.yukawa_type = static_cast<THDM_yukawa_type>(yukawa_type);
   cbasis.mh = basis.mh;
   cbasis.mH = basis.mH;
   cbasis.mA = basis.mA;
   cbasis.mHp = basis.mHp;
   cbasis.sin_beta_minus_alpha = basis.sin_beta_minus_alpha;
   cbasis.lambda_6 = basis.lambda_6;
   cbasis.lambda_7 = basis.lambda_7;
   cbasis.tan_beta = basis.tan_beta;
   cbasis.m122 = basis.m122;
   cbasis.zeta_u = basis.zeta_u;
   cbasis.zeta_d = basis.zeta_d;
   cbasis.zeta_l = basis.zeta_l;
   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         cbasis.Delta_u[i][k] = basis.Delta_u(i, k);
         cbasis.Delta_d[i][k] = basis.Delta_d(i, k);
         cbasis.Delta_l[i][k] = basis.Delta_l(i, k);
         cbasis.Pi_u[i][k] = basis.Pi_u(i, k);
         cbasis.Pi_d[i][k] = basis.Pi_d(i, k);
         cbasis.Pi_l[i][k] = basis.Pi_l(i, k);
      }
   }

   gm2calc::SM cppsm;
   SM csm;
   setup_SM(cppsm, csm);

   gm2calc::thdm::Config cppconfig;
   THDM_config cconfig;
   gm2calc_thdm_config_set_to_default(&cconfig);

   THDM* mc = nullptr;
   gm2calc_thdm_new_with_mass_basis(&mc, &cbasis, &csm, &cconfig);

   return std::make_pair(gm2calc::THDM(basis, cppsm, cppconfig), mc);
}


std::pair<gm2calc::THDM, THDM*> setup_gauge_basis(gm2calc::thdm::Yukawa_type yukawa_type)
{
   gm2calc::thdm::Gauge_basis basis;
   basis.yukawa_type = yukawa_type;
   basis.lambda << 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1;
   basis.tan_beta = 3;
   basis.m122 = 40000;
   basis.zeta_u = 0.1;
   basis.zeta_d = 0.2;
   basis.zeta_l = 0.3;
   basis.Delta_u << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9;
   basis.Delta_d = 2.0*basis.Delta_u;
   basis.Delta_l = 3.0*basis.Delta_u;
   basis.Pi_u << 4.0*basis.Delta_u;
   basis.Pi_d = 2.0*basis.Pi_u;
   basis.Pi_l = 3.0*basis.Pi_u;

   THDM_gauge_basis cbasis;
   cbasis.yukawa_type = static_cast<THDM_yukawa_type>(yukawa_type);
   for (int i = 0; i < 7; i++) {
      cbasis.lambda[i] = basis.lambda(i);
   }
   cbasis.tan_beta = basis.tan_beta;
   cbasis.m122 = basis.m122;
   cbasis.zeta_u = basis.zeta_u;
   cbasis.zeta_d = basis.zeta_d;
   cbasis.zeta_l = basis.zeta_l;
   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         cbasis.Delta_u[i][k] = basis.Delta_u(i, k);
         cbasis.Delta_d[i][k] = basis.Delta_d(i, k);
         cbasis.Delta_l[i][k] = basis.Delta_l(i, k);
         cbasis.Pi_u[i][k] = basis.Pi_u(i, k);
         cbasis.Pi_d[i][k] = basis.Pi_d(i, k);
         cbasis.Pi_l[i][k] = basis.Pi_l(i, k);
      }
   }

   gm2calc::SM cppsm;
   SM csm;
   setup_SM(cppsm, csm);

   gm2calc::thdm::Config cppconfig;
   THDM_config cconfig;
   gm2calc_thdm_config_set_to_default(&cconfig);

   THDM* mc = nullptr;
   gm2calc_thdm_new_with_gauge_basis(&mc, &cbasis, &csm, &cconfig);

   return std::make_pair(gm2calc::THDM(basis, cppsm, cppconfig), mc);
}


void test_mass_basis(gm2calc::thdm::Yukawa_type yukawa_type)
{
   const auto models = setup_mass_basis(yukawa_type);
   const auto mcpp = models.first;
   const auto mc   = models.second;

   CHECK(gm2calc_thdm_calculate_amu_1loop(mc) == gm2calc::calculate_amu_1loop(mcpp));
   CHECK(gm2calc_thdm_calculate_amu_2loop(mc) == gm2calc::calculate_amu_2loop(mcpp));
   CHECK(gm2calc_thdm_calculate_uncertainty_amu_0loop(mc) == gm2calc::calculate_uncertainty_amu_0loop(mcpp));
   CHECK(gm2calc_thdm_calculate_uncertainty_amu_1loop(mc) == gm2calc::calculate_uncertainty_amu_1loop(mcpp));
   CHECK(gm2calc_thdm_calculate_uncertainty_amu_2loop(mc) == gm2calc::calculate_uncertainty_amu_2loop(mcpp));

   gm2calc_thdm_free(mc);
}


void test_gauge_basis(gm2calc::thdm::Yukawa_type yukawa_type)
{
   const auto models = setup_gauge_basis(yukawa_type);
   const auto mcpp = models.first;
   const auto mc   = models.second;

   CHECK(gm2calc_thdm_calculate_amu_1loop(mc) == gm2calc::calculate_amu_1loop(mcpp));
   CHECK(gm2calc_thdm_calculate_amu_2loop(mc) == gm2calc::calculate_amu_2loop(mcpp));
   CHECK(gm2calc_thdm_calculate_uncertainty_amu_0loop(mc) == gm2calc::calculate_uncertainty_amu_0loop(mcpp));
   CHECK(gm2calc_thdm_calculate_uncertainty_amu_1loop(mc) == gm2calc::calculate_uncertainty_amu_1loop(mcpp));
   CHECK(gm2calc_thdm_calculate_uncertainty_amu_2loop(mc) == gm2calc::calculate_uncertainty_amu_2loop(mcpp));

   gm2calc_thdm_free(mc);
}


TEST_CASE("mass_basis")
{
   test_mass_basis(gm2calc::thdm::Yukawa_type::type_1);
   test_mass_basis(gm2calc::thdm::Yukawa_type::type_2);
   test_mass_basis(gm2calc::thdm::Yukawa_type::type_X);
   test_mass_basis(gm2calc::thdm::Yukawa_type::type_Y);
   test_mass_basis(gm2calc::thdm::Yukawa_type::aligned);
   test_mass_basis(gm2calc::thdm::Yukawa_type::general);
}


TEST_CASE("gauge_basis")
{
   test_gauge_basis(gm2calc::thdm::Yukawa_type::type_1);
   test_gauge_basis(gm2calc::thdm::Yukawa_type::type_2);
   test_gauge_basis(gm2calc::thdm::Yukawa_type::type_X);
   test_gauge_basis(gm2calc::thdm::Yukawa_type::type_Y);
   test_gauge_basis(gm2calc::thdm::Yukawa_type::aligned);
   test_gauge_basis(gm2calc::thdm::Yukawa_type::general);
}
