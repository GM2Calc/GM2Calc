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

#include "GeneralTHDM/gm2_2loop_helpers.hpp"
#include "gm2_dilog.hpp"
#include "gm2_ffunctions.hpp"

/**
 * \file gm2_2loop_B.cpp
 *
 * Contains functions necessary to calculate the bosonic THDM
 * contributions for g-2 at the 2-loop level.
 *
 * @todo(alex) Catch divergent case: u = 0
 * @todo(alex) Catch divergent case: u = 1
 * @todo(alex) Catch divergent case: w = 0
 * @todo(alex) Catch divergent case: w = 1/4
 * @todo(alex) Catch divergent case: w = u/4
 */

namespace gm2calc {

namespace general_thdm {

namespace {

const double pi = 3.1415926535897932;
const double pi2 = 9.8696044010893586; // Pi^2

double sqr(double x) noexcept { return x*x; }

double pow3(double x) noexcept { return sqr(x)*x; }

/// Eq.(102), arxiv:1607.06292
double YF1(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto c0 = cw2*(-1 + cw2)*(u + 2*w)/u;

   return
      - 72*c0 - 36*c0*std::log(w)
      + 9*(-8*cw4 - 3*u + 2*cw2*(4 + u))*(u + 2*w)/(2*(u-1)*u)*std::log(u)
      - 9*(3 - 10*cw2 + 8*cw4)*w*(u + 2*w)/((4*w-1)*(u-1))*Phi(w,w,1)
      + 9*(8*cw4 + 3*u - 2*cw2*(4 + u))*w*(u + 2*w)/((4*w-u)*(u-1)*u*u)*Phi(w,w,w)
      ;
}

/// Eq.(121), arxiv:1607.06292
double YFZ(double u, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto u2 = u*u;
   const auto lu = std::log(u);
   const auto li = dilog(1.0 - u);
   const auto phi = Phi(u, 1.0, 1.0);

   const auto z1 = 3*(17 - 48*cw2 + 32*cw4); // Eq.(122)
   const auto z2 = 5 - 12*cw2 + 8*cw4;       // Eq.(123)
   const auto z3 = 3*(1 - 3*cw2 + 2*cw4);    // Eq.(124)

   const double res =
      + z1*u*li
      + z2/(2*u2)*(6*(-4 + u)*u + pi2*(4 + 3*u) + 6*u*(4 + u)*lu
                   - 6*(4 + 3*u)*li + 6*u*(2 + u)*phi)
      + z3*u*(6 + pi2*(-4 + u)*u + 3*lu*(4 + (-4 + u)*u*lu))
      + 12*(-4 + u)*u*li + 6*(-2 + u)*phi;

   return res;
}

/// Eq.(125), arxiv:1607.06292
double YFW(double u, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto cw6 = cw4*cw2;
   const auto u2 = u*u;
   const auto u3 = u2*u;

   const double res =
      - 57.0/2*cw2 - 4*cw6*pi2/u2 + 3*cw4*(32 - 3*pi2)/(4*u)
      + 3*(16*cw6 + 9*cw4*u + 12*cw2*u2 - 19*u3)*dilog(1.0 - u/cw2)/(2*u2)
      + 3*cw2*(16*cw2 + 19*u)*(std::log(cw2/u))/(2*u)
      + 3*(4*cw4 - 50*cw2*u + 19*u2)*Phi(u,cw2,cw2)/(2*(4*cw2-u)*u);

   return res;
}

/// Eq.(105), arxiv:1607.06292
double YF2(double u, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto cw6 = cw4*cw2;
   const auto cw8 = cw4*cw4;
   const auto u2 = u*u;
   const auto u3 = u2*u;

   const double f0 = 3.0/4*cw4*(-640 + 576*cw2 + 7*pi2);     // Eq.(106)
   const double f1 = 96*cw6*(11 - 53*cw2 + 36*cw4);          // Eq.(107)
   const double f2 = -3.0/4*cw2*(-66*cw2 - 48*cw4 + 672*cw6);// Eq.(108)
   const double f3 = -3.0/4*cw2*(109 - 430*cw2 + 120*cw4);   // Eq.(109)
   const double f4 = 96*cw6*(-11 + 9*cw2);                   // Eq.(110)
   const double f5 = 45.0/2*cw4 + 192*cw6;                   // Eq.(111)
   const double f6 = 3.0/4*cw2*(157 + 90*cw2);               // Eq.(112)
   const double f7 = -3.0/4*(18 + 61*cw2);                   // Eq.(113)
   const double f8 = -7 + 61*cw2 - 162*cw4 + 96*cw6;         // Eq.(114)
   const double f9 = 1 - 5*cw2 + 10*cw4;                     // Eq.(115)
   const double f10 = -1728*cw8*(-1 + cw2);                  // Eq.(116)
   const double f11 = 3*cw6*(-899 + 768*cw2);                // Eq.(117)
   const double f12 = 387*cw4 - 363*cw6;                     // Eq.(118)
   const double f13 = 9.0/2*cw2*(57 + 106*cw2);              // Eq.(119)
   const double f14 = -15.0/2*(7 + 45*cw2);                  // Eq.(120)

   const double res =
      + 8*cw6*pi2/u2 + f0/u + 393.0/8*cw2
      + (f1/u + f2 + f3*u)*std::log(cw2)/((4*cw2-1)*(4*cw2-u))
      + (f4/u + f5 + f6*u + f7*u2)*std::log(u)/((u-1)*(4*cw2-u))
      - 3.0/2*(32*cw6/u2 + 21*cw4/u + 15*cw2 - 35*u)*dilog(1.0 - u/cw2)
      + (f8 + f9*u)*9*cw2*(-3 + 4*cw2)/2*Phi(cw2,cw2,1)/(sqr(4*cw2-1)*(u-1))
      + (f10/u2 + f11/u + f12 + f13*u + f14*u2 + 105.0/2*u3)*Phi(u,cw2,cw2)
        /(sqr(4*cw2-u)*(u-1))
      ;

   return YFW(u, cw2) + YFZ(u, cw2) + res;
}

/// Eq.(126), arxiv:1607.06292
double YF3(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto cw6 = cw4*cw2;
   const auto cw8 = cw4*cw4;
   const auto u2 = u*u;
   const auto u3 = u2*u;
   const auto u4 = u2*u2;
   const auto u5 = u4*u;
   const auto w2 = w*w;
   const auto lc = std::log(cw2);
   const auto lu = std::log(u);
   const auto lw = std::log(w);

   // Eq.(127)
   const auto a1 = -9*cw2*u3 + 9*cw2*u2*(3*cw2+w) + 27*cw4*u*(w-cw2)
      + 9*(cw8 - 4*cw6*w + 3*cw4*w2);
   // Eq.(128)
   const auto a2 = 9*cw4*w/2 - 9*u2*(5*cw2 + w) + u*(36*cw4 + 153*cw2*w/4)
      + 9*u3;
   // Eq.(129)
   const auto a3 = 9*cw2*u2 - 9.0/2*cw2*u*(4*cw2 + w);
   // Eq.(130)
   const auto a4 = -9.0/2*u2*w*(2*cw4 + 9*cw2*w + 2*w2)
      + 9.0/8*u*w*(32*cw6 + 13*cw4*w + 35*cw2*w2) + 9*u3*w2;
   // Eq.(131)
   const auto a5 = -9*u3*(cw2 + w) - 9*u*(3*cw6 + 2*cw2*w2)
      + 9*u2*(3*cw4 + 4*cw2*w + w2) + 9.0/2*cw4*(2*cw4 - 6*cw2*w + w2);
   // Eq.(132)
   const auto a6 = -9*u4*(9*cw2 + w) + u*(81*cw6*w - 225*cw8) + 9*cw8*(w - cw2)
      - 9.0/2*u2*(3*cw6 + 37*cw4*w) + u3*(198*cw4 + 72*cw2*w) + 9*u5;
   // Eq.(133)
   const auto a7 = -9*cw2*u4 + 18*cw2*u3*(2*cw2 + w) + 36*u*(cw8 - 2*cw6*w)
      - 9*cw2*u2*(6*cw4 - cw2*w + w2) - 9*cw2*(cw2 - 3*w)*(cw6 - 2*cw4*w + cw2*w2);

   const double res =
      + 9*u*(2*cw2 - u + w)/w
      + (a1*(lu - lc) + 9*cw4*(cw4 - 4*cw2*w + 3*w2)*lc)*(lw - lc)/(2*w2*(cw2-w))
      + a2*lu/(w*(4*cw2-u))
      + a3*lw/(w*(cw2-w))
      + a4*lc/(w2*(4*cw2-u)*(cw2-w))
      + a5/(cw2*w2)*dilog(1.0 - u/cw2)
      + a6/(u*cw2*sqr(4*cw2-u)*(cw2-w))*Phi(u,cw2,cw2)
      + a7/(w2*(cw2-w)*(cw4-2*cw2*(u+w)+sqr(u-w)))*Phi(u,w,cw2)
      ;

   return res;
}

/// Eq.(72), arxiv:1607.06292
double T0(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;

   return
      9/cw4*(u - w)*(cw2*(u - w)*(u + 2*w) - pow3(u - w) + cw4*w)
      /(cw4 + sqr(u - w) - 2*cw2*(u + w))*Phi(u,w,cw2);
}

/// Eq.(73), arxiv:1607.06292
double T1(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;

   return 9/cw4*(u - w)*(cw2*w - sqr(u - w))*dilog(1.0 - u/w);
}

/// Eq.(74), arxiv:1607.06292
double T2(double u, double w, double cw2, double xH, double xA, double xHp, int sgn) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto u2 = u*u;

   const auto f6 = (7 - 14*cw2 + 4*cw4)/(4*cw2*sw2);
   const auto f7 = 1 - 6*cw2 + 4*cw4;
   const auto f8 = (13 - 20*cw2 + 4*cw4)/(cw2*sw2);
   const auto f9 = 7 - 12*cw2 + 8*cw4;

   return std::log(u)*(
      + (6*u2 + cw2*(u - xHp) + 2*cw4*(u - xHp))/(2*(u - w))
      + f6*sqr(u - xHp)*(3*cw4 + 3*cw2*(u - xHp) + sqr(u - xHp))/(cw2*(u - w))
      + sgn*f7*3*u2*(u - xHp)/((xA - xH)*(u - w))
      - f8*3*u*sqr(u - xHp)/(2*(u - w))
      - f9*3*u*(u - xHp)/(2*(u - w))
      );
}

/// Eq.(74), arxiv:1607.06292 with positive sign
double T2p(double u, double w, double cw2, double xH, double xA, double xHp) noexcept
{
   return T2(u, w, cw2, xH, xA, xHp, +1);
}

/// Eq.(74), arxiv:1607.06292 with negative sign
double T2m(double u, double w, double cw2, double xH, double xA, double xHp) noexcept
{
   return T2(u, w, cw2, xH, xA, xHp, -1);
}

/// Eq.(75), arxiv:1607.06292
double T4(double u, double w, double cw2, double xH, double xA) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f5 = cw2*(5 - 16*cw2 + 8*cw4)/sw2;

   return (u - w)*std::log(u)/4*f5*(xA*(3 + 2*xH) - sqr(xA) + 3*xH - sqr(xH) - 3);
}

/// Eq.(76), arxiv:1607.06292
double T5(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f6 = (7 - 14*cw2 + 4*cw4)/(4*cw2*sw2);
   const auto f8 = (13 - 20*cw2 + 4*cw4)/(cw2*sw2);

   return std::log(u)*(
      3.0/2*u + f6/cw2*(pow3(u - w) + 3*cw2*sqr(u - w) + 3*cw4*(u - w))
      - 3.0/2*f8*u*(u - w) - cw2/2 - cw4);
}

/// Eq.(77), arxiv:1607.06292
double T6(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto u2 = u*u;

   return 9.0/2*(
      (u - w)*(u2 - 2*u*w + w*(w - cw2))/cw4*std::log(u/w)*std::log(w/cw2)
      + std::log(cw2)/cw2*(2*u2 + u*(cw2 - 4*w) - w*(cw2 - 2*w)));
}

/// Eq.(78), arxiv:1607.06292
double T7(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f5 = cw2*(5 - 16*cw2 + 8*cw4)/sw2;
   const auto s1 = u + w - 1 + std::sqrt(1 + sqr(u - w) - 2*(u + w));

   return f5*(2*(u + w) - sqr(u - w) - 1)*std::log(s1/(2*std::sqrt(u*w)))
      *(u + w - 1 - 4*u*w/s1);
}

/// Eq.(79), arxiv:1607.06292
double T8(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f6 = (7 - 14*cw2 + 4*cw4)/(4*cw2*sw2);
   const auto s2 = u + w - cw2 + std::sqrt(sqr(u + w - cw2) - 4*u*w);

   return
      2*f6*(4*u*w - sqr(u + w - cw2))*std::log(s2/(2*std::sqrt(u*w)))
      *((u + w)/cw2 - 4*u*w/(cw2*s2) - 1);
}

/// Eq.(103), arxiv:1607.06292
double T9(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto u2 = u*u;
   const auto w2 = w*w;

   return
      - 2*(cw4*w + cw2*(u2 + u*w - 2*w2) - pow3(u-w))*Phi(u,w,cw2)
        /((cw2 - w)*(cw4 - 2*cw2*(u+w) + sqr(u-w)))
      + 2*cw4*(u2 - 4*u*w + 2*w2)*Phi(u,w,w)/(w2*(w-cw2)*(u-4*w))
      - 2*(cw2*u*(u-2*w) + w*sqr(u-w))*dilog(1.0 - u/w)/w2;
}

/// Eq.(104), arxiv:1607.06292
double T10(double u, double w, double cw2) noexcept
{
   const auto u2 = u*u;
   const auto w2 = w*w;
   const auto lwu = std::log(w/u);
   const auto lwc = std::log(w/cw2);

   return
      (u2 - cw2*w - 2*u*w + w2)/(2*(cw2-w))*lwu*lwc
      + cw2*(cw2 + 2*u - 2*w)/(2*(cw2-w))*lwc
      + cw2*u/w*lwu
      + cw2/w*(w-u);
}

/// Eq.(99), arxiv:1607.06292
double fb(double u, double w, double al, double cw2) noexcept
{
   return al*pi/(cw2*(-1.0 + cw2))*(u + 2*w);
}

/// Eq.(100), arxiv:1607.06292
double Fm0(double u, double w, double al, double cw2) noexcept
{
   return 1.0/(al*pi) * YF1(u,w,cw2);
}

/// Eq.(101), arxiv:1607.06292
double Fmp(double u, double w, double al, double cw2) noexcept
{
   return (-9*(-1 + cw2))/(al*pi) * (T9(u,w,cw2)/2 + T10(u,w,cw2));
}

} // anonymous namespace

/**
 * Calculates 2-loop bosonic pure electroweak contributions.
 *
 * Eq (49), arxiv:1607:06292
 */
double amu2L_B_EWadd(const THDM_B_parameters& thdm) noexcept
{
   return 2.3e-11 * thdm.eta * thdm.zetal;
}

/**
 * Calculates 2-loop bosonic non-Yukawa contributions.
 *
 * Eq (71), arxiv:1607:06292
 */
double amu2L_B_nonYuk(const THDM_B_parameters& thdm) noexcept
{
   const auto mw2 = sqr(thdm.mw);
   const auto mz2 = sqr(thdm.mz);
   const auto cw2 = mw2/mz2;
   const auto cw4 = cw2*cw2;
   const auto cw6 = cw4*cw2;
   const auto cw8 = cw4*cw4;
   const auto sw2 = 1.0 - cw2;
   const auto sw4 = sw2*sw2;
   const auto mA2 = sqr(thdm.mA);
   const auto mH2 = sqr(thdm.mh(1));
   const auto mHp2 = sqr(thdm.mHp);
   const auto xH = mH2/mz2;
   const auto xA = mA2/mz2;
   const auto xHp = mHp2/mz2;

   const auto f1 = 7.0/2 - 25.0/(2*cw2) + 4*cw2 - 4*cw4;
   const auto f2 = 2*(17 - 24*cw2 + 56*cw4 - 128*cw6 + 64*cw8);
   const auto f3 = (25 - 32*cw2 + 4*cw4)/(cw2*sw2);
   const auto f4 = 13.0/2 - 15*cw2 + 10*cw4;
   const auto f5 = cw2*(5 - 16*cw2 + 8*cw4)/sw2;

   const double res =
      + (xA - xH)/(xA - xHp)*T2p(xA, xH, cw2, xH, xA, xHp)
      + T2m(xH, xHp, cw2, xH, xA, xHp)
      + (xA - xH)/(xA - xHp)*T4(xA, xHp, cw2, xH, xA)
      + T4(xH, xA, cw2, xH, xA)
      + T5(xHp, xH, cw2)
      + T5(xHp, xA, cw2)
      + T2p(xHp, xH, cw2, xH, xA, xHp)
      + T2p(xHp, xA, cw2, xH, xA, xHp)
      + T6(xA, xHp, cw2)
      + T6(xH, xHp, cw2)
      + T7(xA, xH, cw2)
      + T7(xHp, xHp, cw2)*sqr(1 - 2*cw2)
      + T8(xA, xHp, cw2)
      + T8(xH, xHp, cw2)
      - 16.0/3*cw2*sw2*(1 + 8*cw2 - 8*cw4)
      + 8*cw4*sw4/(5*xHp)
      + f2*xHp
      - f3*sqr(xHp)
      + f1*(sqr(xA) + sqr(xH))
      + f3*xHp*(xA + xH)
      + f4*(xA + xH)
      - f5*xA*xH
      + T1(xA, xHp, cw2)
      + T1(xH, xHp, cw2)
      + T0(xA, xHp, cw2)
      + T0(xH, xHp, cw2)
      ;

   const auto pref = sqr(thdm.alpha/(24*pi*cw2*sw2) * thdm.mm/thdm.mz);

   return pref*res;
}

/**
 * Calculates 2-loop bosonic Yukawa contributions.
 *
 * Eq (52), arxiv:1607:06292
 */
double amu2L_B_Yuk(const THDM_B_parameters& thdm) noexcept
{
   const auto tb = thdm.tb;
   const auto zetal = thdm.zetal;
   const auto lambda5 = thdm.lambda5;
   const auto eta = thdm.eta;
   const auto al = thdm.alpha;
   const auto mhSM2 = sqr(thdm.mhSM);
   const auto mH2 = sqr(thdm.mh(1));
   const auto mHp2 = sqr(thdm.mHp);
   const auto mw2 = sqr(thdm.mw);
   const auto mz2 = sqr(thdm.mz);

   const auto cw2 = mw2/mz2;
   const auto xhSM = mhSM2/mz2;
   const auto xHp = mHp2/mz2;
   const auto xH = mH2/mz2;

   // Eq.(91), arxiv:1607.06292
   const auto a000 = fb(xhSM, xHp, al, cw2)*Fm0(xhSM, xHp, al, cw2);
   // Eq.(92), arxiv:1607.06292
   const auto a0z0 = -fb(xH, 0.0, al, cw2)*(
      Fm0(xH, xHp, al, cw2) + Fmp(xH, xHp, al, cw2));
   // Eq.(93), arxiv:1607.06292
   const auto a500 = Fm0(xhSM, xHp, al, cw2);
   // Eq.(94), arxiv:1607.06292
   const auto a5z0 = -0.5*(
      Fm0(xH, xHp, al, cw2) + Fmp(xH, xHp, al, cw2));
   // Eq.(95), arxiv:1607.06292
   const auto a001 =
      + fb(xH, 0.0, al, cw2)*Fm0(xH, xHp, al, cw2)
      - fb(xhSM, 0.0, al, cw2)*Fm0(xhSM, xHp, al, cw2);
   // Eq.(96), arxiv:1607.06292
   const auto a0z1 =
      - (
         + fb(xH, xHp, al, cw2)*(
            + Fm0(xH, xHp, al, cw2)
            + Fmp(xH, xHp, al, cw2)
         )
         - YF3(xH, xHp, cw2)
         // SM contributions with opposite sign:
         - fb(xhSM, xHp, al, cw2)*(
            + Fm0(xhSM, xHp, al, cw2)
            + Fmp(xhSM, xHp, al, cw2)
         )
         + YF3(xhSM, xHp, cw2)
      )
      + YF2(xH, cw2);
   // Eq.(97), arxiv:1607.06292
   const auto a501 =
      + Fm0(xH, xHp, al, cw2)/2
      - Fm0(xhSM, xHp, al, cw2)/2;
   // Eq.(98), arxiv:1607.06292
   const auto a5z1 =
      - Fm0(xH, xHp, al, cw2)
      - Fmp(xH, xHp, al, cw2)
      + Fm0(xhSM, xHp, al, cw2)
      + Fmp(xhSM, xHp, al, cw2);

   const double res =
      + a000
      + a0z0*(tb - 1.0/tb)*zetal
      + a500*lambda5
      + a5z0*(tb - 1.0/tb)*lambda5*zetal
      + (
         + a001*(tb - 1.0/tb)
         + a0z1*zetal
         + a501*(tb - 1.0/tb)*lambda5
         + a5z1*lambda5*zetal
         )*eta;

   const auto pref = sqr(al/(24*pi*cw2*(1.0 - cw2)) * thdm.mm/thdm.mz);

   return pref*res;
}

/**
 * Calculates the sum of the 2-loop bosonic contributions.
 *
 * Eq (48), arxiv:1607:06292
 */
double amu2L_B(const THDM_B_parameters& thdm) noexcept
{
   return amu2L_B_EWadd(thdm) + amu2L_B_nonYuk(thdm) + amu2L_B_Yuk(thdm);
}

} // namespace general_thdm

} // namespace gm2calc
