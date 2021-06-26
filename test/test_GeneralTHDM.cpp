#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "GeneralTHDM/gm2_2loop_helpers.hpp"

#include <fstream>
#include <iomanip>
#include <limits>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


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
T3 = { -1, -1, 1 }; (* l, d, u *)
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

   CHECK_CLOSE(1e10*amu, -1023.0630847112763, 1e-12);
}


// check data of Fig.8 of arxiv:1607.06292
TEST_CASE("fermionic_figure_8")
{
   const auto sqr = [] (double x) { return x*x; };
   const double pi = 3.1415926535897932;
   const double alpha = 1./137.036;
   const double mm = 0.10565837;
   const double mw2 = sqr(80.379);
   const double mz2 = sqr(91.1876);
   const double cw2 = mw2/mz2;
   const double sw2 = 1 - cw2;
   const double pref = sqr(alpha*mm/(2*pi))/(mw2*sw2);
   Eigen::Matrix<double,3,1> mu2;
   Eigen::Matrix<double,3,1> md2;
   Eigen::Matrix<double,3,1> ml2;
   mu2 << sqr(2.16e-3), sqr(1.270), sqr(173.1);
   md2 << sqr(4.67e-3), sqr(0.093), sqr(4.180);
   ml2 << sqr(511e-6), sqr(mm), sqr(1.7768);

   constexpr int N_steps = 1000;
   const double m_start = 0;
   const double m_stop = 500;

   std::array<double, N_steps + 1> mS{};
   std::array<double, N_steps + 1> amu_H_u{}, amu_A_u{}, amu_Hp_u{}, amu_hH_u{};
   std::array<double, N_steps + 1> amu_H_d{}, amu_A_d{}, amu_Hp_d{}, amu_hH_d{};
   std::array<double, N_steps + 1> amu_H_l{}, amu_A_l{}, amu_Hp_l{}, amu_hH_l{};

   for (int i = 0; i <= N_steps; ++i) {
      mS[i] = m_start + i*(m_stop - m_start)/N_steps;

      for (int g = 0; g < 3; ++g) {
         // Eq.(64), first line, first term for h = H
         amu_H_u[i] += gm2calc::general_thdm::fuS(sqr(mS[i]), mu2(g), mw2, mz2);
         amu_H_d[i] += gm2calc::general_thdm::fdS(sqr(mS[i]), md2(g), mw2, mz2);
         amu_H_l[i] += gm2calc::general_thdm::flS(sqr(mS[i]), ml2(g), mw2, mz2);

         // Eq.(64), first line, first term for h = A
         amu_A_u[i] += gm2calc::general_thdm::fuA(sqr(mS[i]), mu2(g), mw2, mz2);
         amu_A_d[i] += gm2calc::general_thdm::fdA(sqr(mS[i]), md2(g), mw2, mz2);
         amu_A_l[i] += gm2calc::general_thdm::flA(sqr(mS[i]), ml2(g), mw2, mz2);

         // Eq.(64), first line, second term
         amu_Hp_u[i] += gm2calc::general_thdm::fuHp(sqr(mS[i]), md2(g), mu2(g), mw2, mz2);
         amu_Hp_d[i] += gm2calc::general_thdm::fdHp(sqr(mS[i]), md2(g), mu2(g), mw2, mz2);
         amu_Hp_l[i] += gm2calc::general_thdm::flHp(sqr(mS[i]), ml2(g), mw2, mz2);

         // Eq.(64), second line
         amu_hH_u[i] += gm2calc::general_thdm::fuS(sqr(125.0), mu2(g), mw2, mz2)
                      - gm2calc::general_thdm::fuS(sqr(mS[i]), mu2(g), mw2, mz2);
         amu_hH_d[i] += gm2calc::general_thdm::fdS(sqr(125.0), md2(g), mw2, mz2)
                      - gm2calc::general_thdm::fdS(sqr(mS[i]), md2(g), mw2, mz2);
         amu_hH_l[i] += gm2calc::general_thdm::flS(sqr(125.0), ml2(g), mw2, mz2)
                      - gm2calc::general_thdm::flS(sqr(mS[i]), ml2(g), mw2, mz2);
      }
   }

   /*
   std::ofstream ostr("figure_8.txt");

   ostr << "# mS^2"
        << '\t' << "a_u^H" << '\t' << "a_d^H" << '\t' << "a_l^H"
        << '\t' << "a_u^A" << '\t' << "a_d^A" << '\t' << "a_l^A"
        << '\t' << "a_u^Hp" << '\t' << "a_d^Hp" << '\t' << "a_l^Hp"
        << '\t' << "a_u^hH" << '\t' << "a_d^hH" << '\t' << "a_l^hH"
        << '\n';

   for (int i = 0; i <= N_steps; ++i) {
      ostr << std::setprecision(std::numeric_limits<double>::digits10 + 1)
           << mS[i]
           << '\t' << pref*amu_H_u[i] << '\t' << pref*amu_H_d[i] << '\t' << pref*amu_H_l[i]
           << '\t' << pref*amu_A_u[i] << '\t' << pref*amu_A_d[i] << '\t' << pref*amu_A_l[i]
           << '\t' << pref*amu_Hp_u[i] << '\t' << pref*amu_Hp_d[i] << '\t' << pref*amu_Hp_l[i]
           << '\t' << pref*amu_hH_u[i] << '\t' << pref*amu_hH_d[i] << '\t' << pref*amu_hH_l[i]
           << '\n';
   }
   */
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

   // CHECK_CLOSE(1e10*amu, -0.9831098470297721, 1e-12);
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
