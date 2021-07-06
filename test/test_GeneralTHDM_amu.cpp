#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "GeneralTHDM/gm2_1loop_helpers.hpp"
#include "GeneralTHDM/gm2_2loop_helpers.hpp"
#include "read_data.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


namespace {

double sqr(double x) noexcept { return x*x; }

} // anonymous namespace


TEST_CASE("1-loop")
{
   const double zetal = 0.1;
   const double eta = 0.001; // terms of O(eta) are neglected

   gm2calc::general_thdm::THDM_1L_parameters pars;
   pars.alpha = 1./137.036;
   pars.mm = 0.10565837;
   pars.ml << 510.999e-6, pars.mm, 1.7768;
   pars.mw = 80.379;
   pars.mz = 91.1876;
   pars.mhSM = 125.09;
   pars.mh << pars.mhSM, 180.0;
   pars.mA = 250.0;
   pars.mHp = 500.0;
   pars.ylh(1,1) = 1 + eta*zetal;
   pars.ylH(1,1) = -zetal + eta;
   pars.ylA(1,1) = -zetal;
   pars.ylHp(1,1) = pars.ylA(1,1);

   const auto a1L = gm2calc::general_thdm::amu1L(pars);
   const auto a1L_approx = gm2calc::general_thdm::amu1L_approx(pars);

   CHECK_CLOSE(a1L*1e16, a1L_approx*1e16, 1e-10);
}


TEST_CASE("1-loop_approximation")
{
   const double zetal = 0.1;
   const double eta = 0.001; // terms of O(eta) are neglected

   gm2calc::general_thdm::THDM_1L_parameters pars;
   pars.alpha = 1./137.036;
   pars.mm = 0.10565837;
   pars.mw = 80.379;
   pars.mz = 91.1876;
   pars.mhSM = 125.09;
   pars.mh << pars.mhSM, 180.0;
   pars.mA = 250.0;
   pars.mHp = 500.0;
   pars.ylh(1,1) = 1 + eta*zetal;
   pars.ylH(1,1) = -zetal + eta;
   pars.ylA(1,1) = -zetal;
   pars.ylHp(1,1) = pars.ylA(1,1);

   const auto a1L = gm2calc::general_thdm::amu1L_approx(pars);

   const double xH = pars.mh(1)/100.0;
   const double xA = pars.mA/100.0;
   const double xHp = pars.mHp/100.0;

   // Eq.(31), arxiv:1607.06292
   const double a1L_expected =
      sqr(zetal/100.0)*(
         + (3.3 + 0.5*std::log(xH))/sqr(xH)
         - (3.1 + 0.5*std::log(xA))/sqr(xA)
         - 0.04/sqr(xHp)
      );

   CHECK_CLOSE(a1L*1e16, a1L_expected*1e6, 0.01);
}


/* Generated with:

<< ../math/GeneralTHDMBarrZee.m;

MW = 80;
MZ = 90;
cw = MW/MZ;
cw2 = cw^2;
sw2 = 1 - cw2;
MC = 200;
alpha = 1/137;
mmu = 1;
mf1 = {mmu, 2, 3}/10; (* me, md, mu *)
mf2 = {mmu, 2, 3};    (* mm, ms, mc *)
mf3 = {mmu, 2, 3}*10; (* mtau, mb, mt *)
Qf = { -1, -1/3, 2/3 }; (* l, d, u *)
YmuA = 3;
Yf1A = { YmuA, 4, 5 }/10;
Yf2A = { YmuA, 4, 5 };
Yf3A = { YmuA, 4, 5 }*10;

amuC = (
    FCharged[MC, mf1, alpha, mmu, MW, sw2, Qf, Yf1A, YmuA] +
    FCharged[MC, mf2, alpha, mmu, MW, sw2, Qf, Yf2A, YmuA] +
    FCharged[MC, mf3, alpha, mmu, MW, sw2, Qf, Yf3A, YmuA]
)

Print["amuC 10^10 = ", N[10^10 amuC, 17]]

*/

TEST_CASE("2-loop_fermionic_charged")
{
   gm2calc::general_thdm::THDM_F_parameters pars;
   pars.alpha = 1./137;
   pars.mm = 1.0;
   pars.mw = 80.0;
   pars.mz = 90.0;
   pars.mHp = 200.0;
   pars.mu << 0.3, 3.0, 30.0;
   pars.md << 0.2, 2.0, 20.0;
   pars.ml << 0.1, 1.0, 10.0;
   pars.yuS << 0.5, 0.5, 0.5, 5.0, 5.0, 5.0, 50.0, 50.0, 50.0;
   pars.ydS << 0.4, 0.4, 0.4, 4.0, 4.0, 4.0, 40.0, 40.0, 40.0;
   pars.ylS << 0.3, 0.3, 0.3, 3.0, 3.0, 3.0, 30.0, 30.0, 30.0;

   const auto amu = gm2calc::general_thdm::amu2L_F_charged(pars);

   CHECK_CLOSE(1e10*amu, 320.98648556040420, 1e-12);
}

/* Generated with:

MW = 80;
MZ = 90;
cw = MW/MZ;
cw2 = cw^2;
sw2 = 1 - cw2;
MC = 200;
alpha = 1/137;
mmu = 1;
mf1 = {mmu, 2, 3}/10; (* me, md, mu *)
mf2 = {mmu, 2, 3};    (* mm, ms, mc *)
mf3 = {mmu, 2, 3}*10; (* mtau, mb, mt *)
Qf = { -1, -1/3, 2/3 }; (* l, d, u *)
YmuA = 3;
Yf1A = { YmuA, 4, 5 }/10;
Yf2A = { YmuA, 4, 5 };
Yf3A = { YmuA, 4, 5 }*10;
MhSM = 125;
Mphi = { 100, 400, 300 }; (* h, H, A *)
T3 = { -1/2, -1/2, 1/2 }; (* l, d, u *)
gfv = T3/2 - Qf sw2;
Yf1phi = { Yf1A, Yf1A, Yf1A };
Yf2phi = { Yf2A, Yf2A, Yf2A };
Yf3phi = { Yf3A, Yf3A, Yf3A };
Ymuphi = { YmuA, YmuA, YmuA };

amuN = (
    FNeutral[Mphi,mf1,alpha,mmu,MW,MZ,sw2,Qf,gfv,Yf1phi,Ymuphi] -
    Sum[fgammaphi[MhSM,mf1[[j]],alpha,mmu,MW,sw2,Qf[[j]],{1,3,3}[[j]]]
        + fZphi[MhSM,mf1[[j]],alpha,mmu,MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]],
        {j,1,3}] +
    FNeutral[Mphi,mf2,alpha,mmu,MW,MZ,sw2,Qf,gfv,Yf2phi,Ymuphi] -
    Sum[fgammaphi[MhSM,mf2[[j]],alpha,mmu,MW,sw2,Qf[[j]],{1,3,3}[[j]]]
        + fZphi[MhSM,mf2[[j]],alpha,mmu,MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]],
        {j,1,3}] +
    FNeutral[Mphi,mf3,alpha,mmu,MW,MZ,sw2,Qf,gfv,Yf3phi,Ymuphi] -
    Sum[fgammaphi[MhSM,mf3[[j]],alpha,mmu,MW,sw2,Qf[[j]],{1,3,3}[[j]]]
        + fZphi[MhSM,mf3[[j]],alpha,mmu,MW,MZ,sw2,Qf[[j]],{1,3,3}[[j]],gfv[[1]],gfv[[j]]],
        {j,1,3}]
)

Print["amuN 10^10 = ", N[10^10 amuN, 17]]

*/

TEST_CASE("2-loop_fermionic_neutral")
{
   gm2calc::general_thdm::THDM_F_parameters pars;
   pars.alpha = 1./137;
   pars.mm = 1.0;
   pars.mw = 80.0;
   pars.mz = 90.0;
   pars.mh << 100.0, 400.0;
   pars.mhSM = 125.0;
   pars.mA = 300.0;
   pars.mHp = 200.0;
   pars.mu << 0.3, 3.0, 30.0;
   pars.md << 0.2, 2.0, 20.0;
   pars.ml << 0.1, 1.0, 10.0;
   pars.yuS << 0.5, 0.5, 0.5, 5.0, 5.0, 5.0, 50.0, 50.0, 50.0;
   pars.ydS << 0.4, 0.4, 0.4, 4.0, 4.0, 4.0, 40.0, 40.0, 40.0;
   pars.ylS << 0.3, 0.3, 0.3, 3.0, 3.0, 3.0, 30.0, 30.0, 30.0;

   const auto amu = gm2calc::general_thdm::amu2L_F_neutral(pars);

   CHECK_CLOSE(1e10*amu, -683.88944303928634, 1e-12);
}


// check data of Fig.8 of arxiv:1607.06292
TEST_CASE("fermionic_figure_8")
{
   const double pi = 3.1415926535897932;
   const double alpha = 1./137.0;
   const double mm = 0.105658;
   const double mw2 = sqr(80.385);
   const double mz2 = sqr(91.1876);
   const double mh2 = sqr(125.09);
   const double cw2 = mw2/mz2;
   const double sw2 = 1 - cw2;
   const double pref_neut = sqr(alpha*mm/(2*pi))/(mw2*sw2);
   const double pref_char = sqr(alpha*mm/(4*pi))/(2*mw2*sqr(sw2));
   Eigen::Matrix<double,3,1> mu2;
   Eigen::Matrix<double,3,1> md2;
   Eigen::Matrix<double,3,1> ml2;
   mu2 << sqr(2.2e-3), sqr(1.28), sqr(173.34);
   md2 << sqr(4.7e-3), sqr(0.096), sqr(4.18);
   ml2 << sqr(510.999e-6), sqr(mm), sqr(1.7768);

   // Eq.(64), first line, first term for f = u, h = H
   const auto amu_H_u = [&](double mS) {
      return pref_neut*(
         + gm2calc::general_thdm::fuS(sqr(mS), mu2(0), mw2, mz2)
         + gm2calc::general_thdm::fuS(sqr(mS), mu2(1), mw2, mz2)
         + gm2calc::general_thdm::fuS(sqr(mS), mu2(2), mw2, mz2));
   };

   // Eq.(64), first line, first term for f = d, h = H
   const auto amu_H_d = [&](double mS) {
      return pref_neut*(
         + gm2calc::general_thdm::fdS(sqr(mS), md2(0), mw2, mz2)
         + gm2calc::general_thdm::fdS(sqr(mS), md2(1), mw2, mz2)
         + gm2calc::general_thdm::fdS(sqr(mS), md2(2), mw2, mz2));
   };

   // Eq.(64), first line, first term for f = l, h = H
   const auto amu_H_l = [&](double mS) {
      return pref_neut*(
         + gm2calc::general_thdm::flS(sqr(mS), ml2(0), mw2, mz2)
         + gm2calc::general_thdm::flS(sqr(mS), ml2(1), mw2, mz2)
         + gm2calc::general_thdm::flS(sqr(mS), ml2(2), mw2, mz2));
   };

   // Eq.(64), first line, first term for f = u, h = A
   const auto amu_A_u = [&](double mS) {
      return -pref_neut*(
         + gm2calc::general_thdm::fuA(sqr(mS), mu2(0), mw2, mz2)
         + gm2calc::general_thdm::fuA(sqr(mS), mu2(1), mw2, mz2)
         + gm2calc::general_thdm::fuA(sqr(mS), mu2(2), mw2, mz2));
   };

   // Eq.(64), first line, first term for f = d, h = A
   const auto amu_A_d = [&](double mS) {
      return pref_neut*(
         + gm2calc::general_thdm::fdA(sqr(mS), md2(0), mw2, mz2)
         + gm2calc::general_thdm::fdA(sqr(mS), md2(1), mw2, mz2)
         + gm2calc::general_thdm::fdA(sqr(mS), md2(2), mw2, mz2));
   };

   // Eq.(64), first line, first term for f = l, h = A
   const auto amu_A_l = [&](double mS) {
      return pref_neut*(
         + gm2calc::general_thdm::flA(sqr(mS), ml2(0), mw2, mz2)
         + gm2calc::general_thdm::flA(sqr(mS), ml2(1), mw2, mz2)
         + gm2calc::general_thdm::flA(sqr(mS), ml2(2), mw2, mz2));
   };

   // Eq.(64), first line, second term
   const auto amu_Hp_u = [&](double mS) {
      return -pref_char*(
         + gm2calc::general_thdm::fuHp(sqr(mS), md2(0), mu2(0), mw2, mz2)
         + gm2calc::general_thdm::fuHp(sqr(mS), md2(1), mu2(1), mw2, mz2)
         + gm2calc::general_thdm::fuHp(sqr(mS), md2(2), mu2(2), mw2, mz2));
   };

   // Eq.(64), first line, second term
   const auto amu_Hp_d = [&](double mS) {
      return pref_char*(
         + gm2calc::general_thdm::fdHp(sqr(mS), md2(0), mu2(0), mw2, mz2)
         + gm2calc::general_thdm::fdHp(sqr(mS), md2(1), mu2(1), mw2, mz2)
         + gm2calc::general_thdm::fdHp(sqr(mS), md2(2), mu2(2), mw2, mz2));
   };

   // Eq.(64), first line, second term
   const auto amu_Hp_l = [&](double mS) {
      return pref_char*(
         + gm2calc::general_thdm::flHp(sqr(mS), ml2(0), mw2, mz2)
         + gm2calc::general_thdm::flHp(sqr(mS), ml2(1), mw2, mz2)
         + gm2calc::general_thdm::flHp(sqr(mS), ml2(2), mw2, mz2));
   };

   // Eq.(64), second line for f = u
   const auto amu_hH_u = [&](double mS) {
      return pref_neut*(
         + gm2calc::general_thdm::fuS(mh2    , mu2(0), mw2, mz2)
         - gm2calc::general_thdm::fuS(sqr(mS), mu2(0), mw2, mz2)
         + gm2calc::general_thdm::fuS(mh2    , mu2(1), mw2, mz2)
         - gm2calc::general_thdm::fuS(sqr(mS), mu2(1), mw2, mz2)
         + gm2calc::general_thdm::fuS(mh2    , mu2(2), mw2, mz2)
         - gm2calc::general_thdm::fuS(sqr(mS), mu2(2), mw2, mz2)
         );
   };

   // Eq.(64), second line for f = d
   const auto amu_hH_d = [&](double mS) {
      return pref_neut*(
         + gm2calc::general_thdm::fdS(mh2    , md2(0), mw2, mz2)
         - gm2calc::general_thdm::fdS(sqr(mS), md2(0), mw2, mz2)
         + gm2calc::general_thdm::fdS(mh2    , md2(1), mw2, mz2)
         - gm2calc::general_thdm::fdS(sqr(mS), md2(1), mw2, mz2)
         + gm2calc::general_thdm::fdS(mh2    , md2(2), mw2, mz2)
         - gm2calc::general_thdm::fdS(sqr(mS), md2(2), mw2, mz2)
         );
   };

   // Eq.(64), second line, all terms ~ eta*zeta_l
   const auto amu_hH_l = [&](double mS) {
      double res = 0.0;
      for (int g = 0; g < 3; ++g) {
         res += 2 * (gm2calc::general_thdm::flS(mh2, ml2(g), mw2, mz2) -
                     gm2calc::general_thdm::flS(sqr(mS), ml2(g), mw2, mz2)) +
                gm2calc::general_thdm::fuS(mh2    , mu2(g), mw2, mz2) -
                gm2calc::general_thdm::fuS(sqr(mS), mu2(g), mw2, mz2) +
                gm2calc::general_thdm::fdS(mh2    , md2(g), mw2, mz2) -
                gm2calc::general_thdm::fdS(sqr(mS), md2(g), mw2, mz2);
      }
      return pref_neut * res;
   };

   const auto data = gm2calc::test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + PATH_SEPARATOR + "figure_8" +
      PATH_SEPARATOR + "figure_8.txt");

   const double eps = 1e-10;

   for (const auto& p: data) {
      const auto mS = p.at(0);
      const auto a_H_u = p.at(1);
      const auto a_H_d = p.at(2);
      const auto a_H_l = p.at(3);
      const auto a_A_u = p.at(4);
      const auto a_A_d = p.at(5);
      const auto a_A_l = p.at(6);
      const auto a_Hp_u = p.at(7);
      const auto a_Hp_d = p.at(8);
      const auto a_Hp_l = p.at(9);
      const auto a_hH_u = p.at(10);
      const auto a_hH_d = p.at(11);
      const auto a_hH_l = p.at(12);
      CHECK_CLOSE(std::abs(1e10*amu_H_u(mS)), std::abs(1e10*a_H_u), eps);
      CHECK_CLOSE(std::abs(1e10*amu_H_d(mS)), std::abs(1e10*a_H_d), eps);
      CHECK_CLOSE(std::abs(1e10*amu_H_l(mS)), std::abs(1e10*a_H_l), eps);
      CHECK_CLOSE(std::abs(1e10*amu_A_u(mS)), std::abs(1e10*a_A_u), eps);
      CHECK_CLOSE(std::abs(1e10*amu_A_d(mS)), std::abs(1e10*a_A_d), eps);
      CHECK_CLOSE(std::abs(1e10*amu_A_l(mS)), std::abs(1e10*a_A_l), eps);
      CHECK_CLOSE(std::abs(1e10*amu_Hp_u(mS)), std::abs(1e10*a_Hp_u), eps);
      CHECK_CLOSE(std::abs(1e10*amu_Hp_d(mS)), std::abs(1e10*a_Hp_d), eps);
      CHECK_CLOSE(std::abs(1e10*amu_Hp_l(mS)), std::abs(1e10*a_Hp_l), eps);
      CHECK_CLOSE(std::abs(1e10*amu_hH_u(mS)), std::abs(1e10*a_hH_u), eps);
      CHECK_CLOSE(std::abs(1e10*amu_hH_d(mS)), std::abs(1e10*a_hH_d), eps);
      CHECK_CLOSE(std::abs(1e10*amu_hH_l(mS)), std::abs(1e10*a_hH_l), eps);
   }

   /*
   constexpr int N_steps = 1000;
   const double m_start = 0;
   const double m_stop = 500;

   std::ofstream ostr("figure_8.txt");

   // ostr << "# mS^2"
   //      << '\t' << "a_u^H" << '\t' << "a_d^H" << '\t' << "a_l^H"
   //      << '\t' << "a_u^A" << '\t' << "a_d^A" << '\t' << "a_l^A"
   //      << '\t' << "a_u^Hp" << '\t' << "a_d^Hp" << '\t' << "a_l^Hp"
   //      << '\t' << "a_u^hH" << '\t' << "a_d^hH" << '\t' << "a_l^hH"
   //      << '\n';

   for (int i = 0; i < N_steps; ++i) {
      const double mS = m_start + (i + 1.0)*(m_stop - m_start)/N_steps;

      ostr << std::setprecision(std::numeric_limits<double>::digits10 + 1)
           << mS
           << '\t' << amu_H_u(mS) << '\t' << amu_H_d(mS) << '\t' << amu_H_l(mS)
           << '\t' << amu_A_u(mS) << '\t' << amu_A_d(mS) << '\t' << amu_A_l(mS)
           << '\t' << amu_Hp_u(mS) << '\t' << amu_Hp_d(mS) << '\t' << amu_Hp_l(mS)
           << '\t' << amu_hH_u(mS) << '\t' << amu_hH_d(mS) << '\t' << amu_hH_l(mS)
           << '\n';
   }
   */
}


TEST_CASE("2-loop_bosonic_EWadd")
{
   gm2calc::general_thdm::THDM_B_parameters pars;
   pars.alpha = 1/137.036;
   pars.mm = 0.10565837;
   pars.mw = 80.379;
   pars.mz = 91.1876;
   pars.mhSM = 125.0;
   pars.mh << pars.mhSM, 0.0;
   pars.zetal = 2;
   pars.eta = 3;

   const auto amu = gm2calc::general_thdm::amu2L_B_EWadd(pars);

   // Eq.(49), arxiv:1607.06292
   CHECK_CLOSE(1e10*amu, 2.3e-1*pars.eta*pars.zetal, 0.01);
   CHECK_CLOSE(1e10*amu, 2.304734800631431e-1*pars.eta*pars.zetal, 1e-12);
}


TEST_CASE("2-loop_bosonic_nonYuk")
{
   gm2calc::general_thdm::THDM_B_parameters pars;
   // parameter point from Fig.4a, arxiv:1607.06292
   pars.alpha = 1./137.036;
   pars.mm = 0.10565837;
   pars.mw = 80.379;
   pars.mz = 91.1876;
   pars.mh << 0.0, 300.0;
   pars.mA = 10.0;
   pars.mHp = 500.0;

   const auto amu = gm2calc::general_thdm::amu2L_B_nonYuk(pars);

   CHECK_CLOSE(1e10*amu, -0.9831098470297721, 1e-9);
}


TEST_CASE("2-loop_bosonic_Yuk")
{
   gm2calc::general_thdm::THDM_B_parameters pars;
   pars.alpha = 1./137.036;
   pars.mm = 0.10565837;
   pars.mw = 80.379;
   pars.mz = 91.1876;
   pars.mh << 0.0, 300.0;
   pars.mhSM = 125.0;
   pars.mHp = 500.0;
   pars.tb = 30.0;
   pars.zetal = 1.0;
   pars.eta = 2.0;
   pars.lambda5 = 3.0;

   const auto amu = gm2calc::general_thdm::amu2L_B_Yuk(pars);

   /* number to compare to has been generated with:

pp = {
    AL -> 1/137.036,
    MM -> 0.10565837,
    MW -> 80.379,
    MZ -> 91.1876,
    MH -> 125,
    MHH -> 300,
    MA0 -> 1,
    MHp -> 500,
    CW -> MW/MZ,
    TB -> 30,
    ZetaL -> 1,
    aeps -> 2, (* eta *)
    Lambda5 -> 3
}

10^10 aYuk /. TF[x___] :> TFex[x] //. pp

    */
   CHECK_CLOSE(1e10*amu, -0.5388969587184598, 1e-12);
}
