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
#include "gm2_numerics.hpp"

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

/// Eq.(102), arxiv:1607.06292
double YF1(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto c0 = cw2*(-1 + cw2)*(u + 2*w)/u;

   return
      - 72*c0 - 36*c0*std::log(w)
      + 9*(-8*cw4 - 3*u + 2*cw2*(4 + u))*(u + 2*w)/(2*(u-1)*u)*std::log(u)
      - 9*(3 - 10*cw2 + 8*cw4)*w*(u + 2*w)/((4*w-1)*(u-1))*Phi(w,w,1)
      + 9*(8*cw4 + 3*u - 2*cw2*(4 + u))*w*(u + 2*w)/((4*w-u)*(u-1)*u*u)*Phi(u,w,w)
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
      + z3*u*(6 + pi2*(-4 + u)*u + 3*lu*(4 + (-4 + u)*u*lu)
              + 12*(-4 + u)*u*li + 6*(-2 + u)*phi);

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
      9/cw4*(u - w)*(cw2*(u - w)*(u + 2*w) - cube(u - w) + cw4*w)
      /(cw4 + sqr(u - w) - 2*cw2*(u + w))*Phi(u,w,cw2);
}

/// Eq.(73), arxiv:1607.06292
double T1(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;

   return 9/cw4*(u - w)*(cw2*w - sqr(u - w))*dilog(1.0 - u/w);
}

/// calculates potentially divergent (a^2 Log[a] - b^2 Log[b])/(a - b)
double dxlog(double a, double b) noexcept
{
   const double lb = std::log(b);

   if (is_equal_rel(a, b, 1e-4)) {
      return b*(1 + 2*lb) + (a - b)*(1.5 + lb) + sqr(a - b)/(3*b);
   }

   return (sqr(a)*std::log(a) - sqr(b)*lb)/(a - b);
}

/**
 * Calculates the following combination of T2 loop functions, Eq.(74),
 * arxiv:1607.06292:
 *
 * (xA - xH)/(xA - xHp)*T2p[xA, xH] + T2m[xH, xHp] + T2p[xHp, xH] + T2p[xHp, xA]
 */
double TX(double xH, double xA, double xHp, double cw2) noexcept
{
   const auto cw4 = sqr(cw2);
   const auto sw2 = 1.0 - cw2;
   const auto f6 = (7 - 14*cw2 + 4*cw4)/(4*cw2*sw2);
   const auto f7 = 1 - 6*cw2 + 4*cw4;
   const auto f8 = (13 - 20*cw2 + 4*cw4)/(cw2*sw2);
   const auto f9 = 7 - 12*cw2 + 8*cw4;
   const auto lH = std::log(xH);
   const auto lA = std::log(xA);
   const auto xAH  = dxlog(xA, xH);
   const auto xAHp = dxlog(xA, xHp);
   const auto xHHp = dxlog(xH, xHp);

   return
      (cw2 + 2*cw4 - 3*f9*xA)*lA/2
      + (3*cw2*f6 - (3*f8*xA)/2)*(xA - xHp)*lA
      + 3*f6*sqr(xA - xHp)*lA
      + (f6*cube(xA - xHp)*lA)/cw2
      + (cw2 + 2*cw4 - 3*f9*xH)*lH/2
      + (3*cw2*f6 - (3*f8*xH)/2)*(xH - xHp)*lH
      + 3*f6*sqr(xH - xHp)*lH
      + (f6*cube(xH - xHp)*lH)/cw2
      + 3*(f7*xAH + xAHp + xHHp);
}

/// Eq.(75), arxiv:1607.06292, prefactor (u-v) has been pulled out
double T4(double u, double cw2, double xH, double xA) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f5 = cw2*(5 - 16*cw2 + 8*cw4)/sw2;

   return std::log(u)/4*f5*(xA*(3 + 2*xH) - sqr(xA) + 3*xH - sqr(xH) - 3);
}

/// Eq.(76), arxiv:1607.06292
double T5(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f6 = (7 - 14*cw2 + 4*cw4)/(4*cw2*sw2);
   const auto f8 = (13 - 20*cw2 + 4*cw4)/(cw2*sw2);

   return std::log(u)*(
      3.0/2*u + f6/cw2*(cube(u - w) + 3*cw2*sqr(u - w) + 3*cw4*(u - w))
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

/**
 * Eq (78), arxiv:1607.06292
 *
 * @note overall factor -1/2 is missing in arxiv:1607.06292v2
 */
double T7(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f5 = cw2*(5 - 16*cw2 + 8*cw4)/sw2;
   const auto ra = std::complex<double>(1 + sqr(u - w) - 2*(u + w), 0.0);
   const auto s1 = u + w - 1.0 + std::sqrt(ra); // Eq.(79)

   const auto res =
      -0.5*f5*(2*(u + w) - sqr(u - w) - 1)*gm2calc::log(s1/(2*std::sqrt(u*w)))
      *(u + w - 1 - 4*u*w/s1);

   return std::real(res);
}

/// Eq.(80), arxiv:1607.06292
double T8(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto f6 = (7 - 14*cw2 + 4*cw4)/(4*cw2*sw2);
   const auto ra = std::complex<double>(sqr(u + w - cw2) - 4*u*w, 0.0);
   const auto s2 = u + w - cw2 + std::sqrt(ra); // Eq.(81)

   const auto res =
      2.0*f6*(4*u*w - sqr(u + w - cw2))*gm2calc::log(s2/(2*std::sqrt(u*w)))
      *((u + w)/cw2 - 4.0*u*w/(cw2*s2) - 1.0);

   return std::real(res);
}

/// Eq.(103), arxiv:1607.06292
double T9(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;
   const auto u2 = u*u;
   const auto w2 = w*w;

   return
      - 2*(cw4*w + cw2*(u2 + u*w - 2*w2) - cube(u-w))*Phi(u,w,cw2)
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
   return 1.0/(al*pi) * cw2*(-1 + cw2)/(u + 2*w) * YF1(u,w,cw2);
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
   const auto mw2 = sqr(thdm.mw);
   const auto mz2 = sqr(thdm.mz);
   const auto cw2 = mw2/mz2;
   const auto cw4  = cw2*cw2;
   const auto cw6  = cw4*cw2;
   const auto cw8  = cw4*cw4;
   const auto cw10 = cw8*cw2;
   const auto cw12 = cw8*cw4;
   const auto cw14 = cw8*cw6;
   const auto mh2 = sqr(thdm.mh(0));
   const auto xh = mh2/mz2;
   const auto s0 = std::sqrt(std::complex<double>(1 - 4*cw2, 0.0));
   const auto s1 = std::sqrt(std::complex<double>(xh*(-4 + xh), 0.0));
   const auto s2 = std::sqrt(std::complex<double>(xh*(-4*cw2 + xh), 0.0));
   const auto lh = std::log(xh);
   const auto lc = std::log(cw2);
   const auto l1 = gm2calc::log(1.0 - (s1 + xh)/2.0);
   const auto l2 = gm2calc::log((-s1 + xh)/2.0);
   const auto l3 = gm2calc::log((1.0 - s0)/2.0);
   const auto l4 = gm2calc::log(1.0 - (s2 + xh)/(2.0*cw2));
   const auto l5 = gm2calc::log((-s2 + xh)/(2.0*cw2));
   const auto li1 = dilog((xh - s1)/2.0);
   const auto li2 = dilog(1.0 - (xh + s1)/2.0);
   const auto li3 = dilog(1.0 - xh/cw2);
   const auto li4 = dilog(1.0 - xh);
   const auto li5 = dilog((1.0 - s0)/2.0);
   const auto li6 = dilog(1.0 - (xh + s2)/(2.0*cw2));
   const auto li7 = dilog((xh - s2)/(2.0*cw2));
   const auto l32 = l3*l3;
   const auto lc2 = lc*lc;
   const auto lh2 = lh*lh;

   const auto xm2 =
      -1280.*cw4*(-6.*li4 + pi2) + 8192.*cw6*(-6.*li4 + pi2) -
       768.*cw10*(-4.*li3 + 64.*li4 - 10.*pi2 + 90.*l4*l5*s2 - 90.*li6*s2 - 90.*li7*s2 + 15.*pi2*s2) +
       256.*cw8*(336.*li4 + 54.*(l4*l5 - li6 - li7)*s2 + pi2*(-56. + 9.*s2)) +
       1024.*cw12*(-12.*li3 + 54.*(l4*l5 - li6 - li7)*s2 + pi2*(2. + 9.*s2));

   const auto xm1 =
      2048.*cw12*(108. + 54.*lc - 54.*lh + 6.*li3 - pi2) + 640.*cw2*(-6.*li4 + pi2) -
       64.*cw4*(-120. + 120.*lh - 354.*li4 + 59.*pi2 + 60.*l1*l2*s1 - 60.*li1*s1 - 60.*li2*s1 + 10.*pi2*s1) -
       16.*cw8*(-9024. - 1920.*lc + 7296.*lh - 48.*li3 - 192.*li4 + 40.*pi2 + 2688.*l1*l2*s1 - 2688.*li1*s1 -
          2688.*li2*s1 + 448.*pi2*s1 - 6594.*l4*l5*s2 + 6594.*li6*s2 + 6594.*li7*s2 - 1099.*pi2*s2) +
       4.*cw6*(-12288. + 12288.*lh - 7680.*li4 + 1280.*pi2 + 6144.*l1*l2*s1 - 6144.*li1*s1 - 6144.*li2*s1 +
          1024.*pi2*s1 - 5442.*l4*l5*s2 + 5442.*li6*s2 + 5442.*li7*s2 - 907.*pi2*s2) +
       1024.*cw10*(-330. - 147.*lc + 195.*lh - 6.*li3 + 12.*li4 - pi2 + 24.*l1*l2*s1 - 24.*li1*s1 -
          24.*li2*s1 + 4.*pi2*s1 - 72.*l4*l5*s2 + 72.*li6*s2 + 72.*li7*s2 - 12.*pi2*s2);

   const auto x0 =
    -80.*(-6.*li4 + pi2) + (36864.*cw14*(6.*l32 - 3.*lc2 - 12.*li5 + pi2))/s0 -
    (6912.*cw12*(78.*l32 - 39.*lc2 - 156.*li5 + 13.*pi2 + 32.*s0 + 16.*lc*s0))/s0 +
    32.*cw2*(-120. + 120.*lh - 66.*li4 + 11.*pi2 + 60.*l1*l2*s1 - 60.*li1*s1 - 60.*li2*s1 + 10.*pi2*s1) -
    (64.*cw10*(-6570.*l32 + 3285.*lc2 + 13140.*li5 - 1095.*pi2 - 3744.*s0 - 1668.*lc*s0 + 336.*lh*s0 -
         48.*li3*s0 - 576.*li4*s0 + 104.*pi2*s0 + 192.*l1*l2*s0*s1 - 192.*li1*s0*s1 - 192.*li2*s0*s1 +
         32.*pi2*s0*s1))/s0 - 4.*cw4*
     (-3744. + 4704.*lh + 1632.*li4 - 272.*pi2 + 2592.*l1*l2*s1 - 2592.*li1*s1 - 2592.*li2*s1 + 432.*pi2*s1 -
       1386.*l4*l5*s2 + 1386.*li6*s2 + 1386.*li7*s2 - 231.*pi2*s2) +
    (4.*cw6*(3024.*l32 - 1512.*lc2 - 6048.*li5 + 504.*pi2 - 3348.*s0 - 2244.*lc*s0 + 4356.*lh*s0 -
         96.*li3*s0 + 11136.*li4*s0 - 1816.*pi2*s0 + 2304.*l1*l2*s0*s1 - 2304.*li1*s0*s1 - 2304.*li2*s0*s1 +
         384.*pi2*s0*s1 - 6222.*l4*l5*s0*s2 + 6222.*li6*s0*s2 + 6222.*li7*s0*s2 - 1037.*pi2*s0*s2))/
     s0 + (16.*cw8*(-7596.*l32 + 3798.*lc2 + 15192.*li5 - 1266.*pi2 + 852.*s0 + 1224.*lc*s0 -
         564.*lh*s0 + 48.*li3*s0 - 4416.*li4*s0 + 704.*pi2*s0 + 576.*l1*l2*s0*s1 - 576.*li1*s0*s1 -
         576.*li2*s0*s1 + 96.*pi2*s0*s1 + 678.*l4*l5*s0*s2 - 678.*li6*s0*s2 - 678.*li7*s0*s2 + 113.*pi2*s0*s2))
      /s0;

   const auto x1 =
      (-14592.*cw12*(6.*l32 - 3.*lc2 - 12.*li5 + pi2))/s0 -
       20.*(-24. + 24.*lh + 6.*li4 - pi2 + 12.*l1*l2*s1 - 12.*li1*s1 - 12.*li2*s1 + 2.*pi2*s1) -
       (64.*cw10*(-3762.*l32 + 1881.*lc2 + 7524.*li5 - 627.*pi2 - 1824.*s0 - 684.*lc*s0 - 384.*lh*s0 -
            768.*li4*s0 + 768.*l1*l2*s0*s1 - 768.*li1*s0*s1 - 768.*li2*s0*s1 + 128.*pi2*s0*s1))/s0 +
       (64.*cw8*(-3114.*l32 + 1557.*lc2 + 6228.*li5 - 519.*pi2 - 2853.*s0 - 768.*lc*s0 - 291.*lh*s0 -
            24.*li3*s0 - 1632.*li4*s0 + 58.*pi2*s0 + 1440.*l1*l2*s0*s1 - 1440.*li1*s0*s1 - 1440.*li2*s0*s1 +
            240.*pi2*s0*s1))/s0 + 2.*cw2*
        (864. + 96.*lh + 1824.*li4 - 304.*pi2 + 288.*l1*l2*s1 - 288.*li1*s1 - 288.*li2*s1 + 48.*pi2*s1 +
          270.*l4*l5*s2 - 270.*li6*s2 - 270.*li7*s2 + 45.*pi2*s2) +
       (4.*cw4*(-1512.*l32 + 756.*lc2 + 3024.*li5 - 252.*pi2 - 5190.*s0 - 345.*lc*s0 + 3081.*lh*s0 -
            804.*li3*s0 - 6576.*li4*s0 + 818.*pi2*s0 + 2496.*l1*l2*s0*s1 - 2496.*li1*s0*s1 - 2496.*li2*s0*s1 +
            416.*pi2*s0*s1 - 198.*l4*l5*s0*s2 + 198.*li6*s0*s2 + 198.*li7*s0*s2 - 33.*pi2*s0*s2))/s0
        - (16.*cw6*(-3690.*l32 + 1845.*lc2 + 7380.*li5 - 615.*pi2 - 3939.*s0 - 780.*lc*s0 + 1158.*lh*s0 -
            828.*li3*s0 - 4848.*li4*s0 + 348.*pi2*s0 + 3360.*l1*l2*s0*s1 - 3360.*li1*s0*s1 - 3360.*li2*s0*s1 +
            560.*pi2*s0*s1 + 342.*l4*l5*s0*s2 - 342.*li6*s0*s2 - 342.*li7*s0*s2 + 57.*pi2*s0*s2))/s0;

   const auto x2 =
      (384.*cw10*(6.*l32 - 3.*lc2 - 12.*li5 + pi2 - 48.*s0 - 96.*lh*s0 - 96.*lh2*s0 - 512.*li4*s0 -
            32.*pi2*s0 + 144.*l1*l2*s0*s1 - 144.*li1*s0*s1 - 144.*li2*s0*s1 + 24.*pi2*s0*s1))/s0 +
       (12.*cw6*(1734.*l32 - 867.*lc2 - 3468.*li5 + 289.*pi2 + 1368.*s0 - 232.*lc*s0 - 964.*lh*s0 -
            2688.*lh2*s0 - 1072.*li3*s0 - 10688.*li4*s0 - 936.*pi2*s0 + 384.*l1*l2*s0*s1 - 384.*li1*s0*s1 -
            384.*li2*s0*s1 + 64.*pi2*s0*s1))/s0 -
       (16.*cw8*(1206.*l32 - 603.*lc2 - 2412.*li5 + 201.*pi2 - 960.*s0 + 72.*lc*s0 - 3264.*lh*s0 -
            4032.*lh2*s0 - 19968.*li4*s0 - 1344.*pi2*s0 + 4512.*l1*l2*s0*s1 - 4512.*li1*s0*s1 -
            4512.*li2*s0*s1 + 752.*pi2*s0*s1))/s0 +
       cw2*(3867. + 426.*lc - 2106.*lh + 1572.*li3 + 5568.*li4 - 384.*pi2 + (756.*l32)/s0 -
          (378.*lc2)/s0 - (1512.*li5)/s0 + (126.*pi2)/s0 - 4032.*l1*l2*s1 +
          4032.*li1*s1 + 4032.*li2*s1 - 672.*pi2*s1 - 420.*l4*l5*s2 + 420.*li6*s2 + 420.*li7*s2 - 70.*pi2*s2)
        + 4.*(-150. + 90.*lh - 90.*li4 + 15.*pi2 + 30.*l1*l2*s1 - 30.*li1*s1 - 30.*li2*s1 + 5.*pi2*s1 -
          48.*l4*l5*s2 + 48.*li6*s2 + 48.*li7*s2 - 8.*pi2*s2) +
       (6.*cw4*(-1122.*l32 + 561.*lc2 + 2244.*li5 - 187.*pi2 - 1774.*s0 - 48.*lc*s0 + 478.*lh*s0 +
            768.*lh2*s0 - 512.*li3*s0 - 224.*li4*s0 + 372.*pi2*s0 + 2784.*l1*l2*s0*s1 - 2784.*li1*s0*s1 -
            2784.*li2*s0*s1 + 464.*pi2*s0*s1 + 792.*l4*l5*s0*s2 - 792.*li6*s0*s2 - 792.*li7*s0*s2 +
            132.*pi2*s0*s2))/s0;

   const auto x3 =
      -3072.*cw10*(-15.*lh2 - 60.*li4 - 5.*pi2 + 6.*l1*l2*s1 - 6.*li1*s1 - 6.*li2*s1 + pi2*s1) +
       (2.*cw4*(342.*l32 - 171.*lc2 - 684.*li5 + 57.*pi2 + 3366.*s0 + 834.*lc*s0 + 6444.*lh*s0 +
            5184.*lh2*s0 + 3144.*li3*s0 + 29184.*li4*s0 + 1728.*pi2*s0 - 8256.*l1*l2*s0*s1 + 8256.*li1*s0*s1 +
            8256.*li2*s0*s1 - 1376.*pi2*s0*s1))/s0 +
       (48.*cw8*(30.*l32 - 15.*lc2 - 60.*li5 + 5.*pi2 + 192.*s0 + 384.*lh*s0 - 1296.*lh2*s0 - 4672.*li4*s0 -
            432.*pi2*s0 + 96.*l1*l2*s0*s1 - 96.*li1*s0*s1 - 96.*li2*s0*s1 + 16.*pi2*s0*s1))/s0 +
       (4.*cw6*(-450.*l32 + 225.*lc2 + 900.*li5 - 75.*pi2 - 3936.*s0 - 180.*lc*s0 - 7680.*lh*s0 +
            2016.*lh2*s0 - 1920.*li4*s0 + 672.*pi2*s0 + 7296.*l1*l2*s0*s1 - 7296.*li1*s0*s1 -
            7296.*li2*s0*s1 + 1216.*pi2*s0*s1))/s0 +
       4.*(-6. - 15.*lh - 48.*li3 - 102.*li4 + 102.*l1*l2*s1 - 102.*li1*s1 - 102.*li2*s1 + 17.*pi2*s1 +
          48.*l4*l5*s2 - 48.*li6*s2 - 48.*li7*s2 + 8.*pi2*s2) -
       (1.*cw2*(108.*l32 - 54.*lc2 - 216.*li5 + 18.*pi2 + 747.*s0 + 426.*lc*s0 + 1350.*lh*s0 + 2304.*lh2*s0 +
            804.*li3*s0 + 9696.*li4*s0 + 768.*pi2*s0 - 672.*l1*l2*s0*s1 + 672.*li1*s0*s1 + 672.*li2*s0*s1 -
            112.*pi2*s0*s1 + 768.*l4*l5*s0*s2 - 768.*li6*s0*s2 - 768.*li7*s0*s2 + 128.*pi2*s0*s2))/s0;

   const auto x4 =
      -3072.*cw10*(3.*lh2 + 12.*li4 + pi2) +
       768.*cw8*(-9.*lh2 - 36.*li4 - 3.*pi2 + 12.*l1*l2*s1 - 12.*li1*s1 - 12.*li2*s1 + 2.*pi2*s1) -
       24.*(-6. - 12.*lh - 12.*lh2 - 8.*li3 - 65.*li4 - 4.*pi2 + 18.*l1*l2*s1 - 18.*li1*s1 - 18.*li2*s1 +
          3.*pi2*s1) + 48.*cw4*(42. + 84.*lh - 312.*lh2 - 1136.*li4 - 104.*pi2 + 42.*l1*l2*s1 - 42.*li1*s1 -
          42.*li2*s1 + 7.*pi2*s1) - 192.*cw6*
        (6. + 12.*lh - 156.*lh2 - 608.*li4 - 52.*pi2 + 66.*l1*l2*s1 - 66.*li1*s1 - 66.*li2*s1 + 11.*pi2*s1) +
       24.*cw2*(-42. - 84.*lh + 36.*lh2 - 32.*li3 + 28.*li4 + 12.*pi2 + 78.*l1*l2*s1 - 78.*li1*s1 -
          78.*li2*s1 + 13.*pi2*s1);

   const auto x5 =
      1536.*cw8*(3.*lh2 + 12.*li4 + pi2) +
       24.*(-15.*lh2 - 60.*li4 - 5.*pi2 + 6.*l1*l2*s1 - 6.*li1*s1 - 6.*li2*s1 + pi2*s1) +
       336.*cw4*(-3.*lh2 - 12.*li4 - pi2 + 6.*l1*l2*s1 - 6.*li1*s1 - 6.*li2*s1 + pi2*s1) -
       192.*cw6*(27.*lh2 + 108.*li4 + 9.*pi2 + 6.*l1*l2*s1 - 6.*li1*s1 - 6.*li2*s1 + pi2*s1) -
       24.*cw2*(-81.*lh2 - 324.*li4 - 27.*pi2 + 42.*l1*l2*s1 - 42.*li1*s1 - 42.*li2*s1 + 7.*pi2*s1);

   const auto x6 =
      24.*(3.*lh2 + 12.*li4 + pi2) - 168.*cw2*(3.*lh2 + 12.*li4 + pi2) + 336.*cw4*(3.*lh2 + 12.*li4 + pi2) -
       192.*cw6*(3.*lh2 + 12.*li4 + pi2);

   const auto res = (xm2/xh + xm1)/xh + x0 + xh*(x1 + xh*(x2 + xh*(x3 + xh*(x4 + xh*(x5 + xh*x6)))));

   const auto pref = sqr(thdm.alpha_em*thdm.mm/(48*thdm.mz*pi*cw2*(4*cw2 - xh)*(1 - cw2)))
      /(2*(-1 + 4*cw2)*(-1 + xh));

   return pref * std::real(res) * thdm.eta * thdm.zetal;
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
      + TX(xH, xA, xHp, cw2)
      + (xA - xH)*(T4(xA, cw2, xH, xA) - T4(xH, cw2, xH, xA))
      + T5(xHp, xH, cw2)
      + T5(xHp, xA, cw2)
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

   const auto pref = sqr(thdm.alpha_em/(24*pi*cw2*sw2) * thdm.mm/thdm.mz);

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
   const auto al = thdm.alpha_em;
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
