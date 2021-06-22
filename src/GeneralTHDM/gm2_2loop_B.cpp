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
 */

namespace gm2calc {

namespace general_thdm {

namespace {

const double pi = 3.1415926535897932;

double sqr(double x) noexcept { return x*x; }

double pow3(double x) noexcept { return sqr(x)*x; }

/// Eq.(102), arxiv:1607.06292
double YF1(double u, double w, double cw2) noexcept
{
   const auto cw4 = cw2*cw2;

   // @todo(alex) avoid re-calculation of common sub-expressions
   return
      - 72*cw2*(-1 + cw2)*(u + 2*w)/u
      - 36*cw2*(-1 + cw2)*(u + 2*w)/u*std::log(w)
      + 9*(-8*cw4 - 3*u + 2*cw2*(4 + u))*(u + 2*w)/(2*(u-1)*u)*std::log(u)
      - 9*(3 - 10*cw2 + 8*cw4)*w*(u + 2*w)/((4*w-1)*(u-1))*Phi(w,w,1)
      + 9*(8*cw4 + 3*u - 2*cw2*(4 + u))*w*(u + 2*w)/((4*w-u)*(u-1)*u*u)*Phi(w,w,w)
      ;
}

/// Eq.(121), arxiv:1607.06292
double YFZ(double u, double w, double al, double cw2, double mm2, double mz2) noexcept
{
   const auto al2 = al*al;
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto sw4 = sw2*sw2;
   const auto u2 = u*u;
   const auto u4 = u2*u2;

   const auto z1 = 3*(17 - 48*cw2 + 32*cw4); // Eq.(122)
   const auto z2 = 5 - 12*cw2 + 8*cw4;       // Eq.(123)
   const auto z3 = 3*(1 - 3*cw2 + 2*cw4);    // Eq.(124)

   // @todo(alex) avoid re-calculation of common sub-expressions
   const double res =
      + z1*u*dilog(1.0 - u)
      + z2/(2*u2)*(6*(-4 + u)*u + sqr(pi)*(4 + 3*u) + 6*u*(4 + u)*std::log(u)
                   - 6*(4 + 3*u)*dilog(1.0 - u) + 6*u*(2 + u)*Phi(u,1.0,1.0))
      + z3*u*(6 + sqr(pi)*(-4 + u)*u + 3*std::log(u)*(4 + (-4 + u)*u*std::log(u)))
      + 12*(-4 + u)*u*dilog(1.0 - u) + 6*(-2 + u)*Phi(u,1.0,1.0)
      ;

   return al2/(576*sqr(pi)*cw4*sw4) * mm2/mz2 * res;
}

/// Eq.(126), arxiv:1607.06292
double YF3(double u, double w, double al, double cw2, double mm2, double mz2) noexcept
{
   const auto al2 = al*al;
   const auto cw4 = cw2*cw2;
   const auto cw6 = cw4*cw2;
   const auto cw8 = cw4*cw4;
   const auto sw2 = 1.0 - cw2;
   const auto sw4 = sw2*sw2;
   const auto u2 = u*u;
   const auto u3 = u2*u;
   const auto u4 = u2*u2;
   const auto u5 = u4*u;
   const auto w2 = w*w;

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

   // @todo(alex) avoid re-calculation of common sub-expressions
   const double res =
      + 9*u*(2*cw2 - u + w)/w
      + (a1*std::log(u/cw2) + 9*cw4*(cw4 - 4*cw2*w + 3*w2)*std::log(cw2))
        *std::log(w/cw2)/(2*w2*(cw2-w))
      + a2*std::log(u)/(w*(4*cw2-u))
      + a3*std::log(w)/(w*(cw2-w))
      + a4*std::log(cw2)/(w2*(4*cw2-u)*(cw2-w))
      + a5/(cw2*w2)*dilog(1.0 - u/cw2)
      + a6/(u*cw2*sqr(4*cw2-u)*(cw2-w))*Phi(u,cw2,cw2)
      + a7/(w2*(cw2-w)*(cw4-2*cw2*(u+w)+sqr(u-w)))*Phi(u,w,cw2)
      ;

   return al2/(576*sqr(pi)*cw4*sw4) * mm2/mz2 * res;
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

   // @todo(alex) avoid re-calculation of common sub-expressions
   return
      (u2 - cw2*w - 2*u*w + w2)/(2*(cw2-w))*std::log(w/u)*std::log(w/cw2)
      + cw2*(cw2 + 2*u - 2*w)/(2*(cw2-w))*std::log(w/cw2)
      + cw2*u/w*std::log(w/u)
      + cw2/w*(w-u);
}

/// Eq.(99), arxiv:1607.06292
double fb(double u, double w, double al, double cw2) noexcept
{
   return al*pi/(cw2*(-1.0 + cw2))*(u + 2*w);
}

/// Eq.(100), arxiv:1607.06292
double Fm0(double u, double w, double al, double cw2, double mm2, double mz2) noexcept
{
   const auto al2 = al*al;
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto sw4 = sw2*sw2;

   // @todo(alex) cancel out al*pi
   return al2/(576*sqr(pi)*cw4*sw4) * mm2/mz2
      * 1.0/(al*pi) * YF1(u,w,cw2);
}

/// Eq.(101), arxiv:1607.06292
double Fmp(double u, double w, double al, double cw2, double mm2, double mz2) noexcept
{
   const auto al2 = al*al;
   const auto cw4 = cw2*cw2;
   const auto sw2 = 1.0 - cw2;
   const auto sw4 = sw2*sw2;

   // @todo(alex) cancel out al*pi
   return al2/(576*sqr(pi)*cw4*sw4) * mm2/mz2
      * (-9*(-1 + cw2))/(al*pi) * (T9(u,w,cw2)/2 + T10(u,w,cw2));
}

} // anonymous namespace

/**
 * Calculates 2-loop bosonic pure electroweak contributions.
 *
 * Eq (49), arxiv:1607:06292
 */
double amu2L_B_EWadd(double eta, double zetal)
{
   return 2.3e-11 * eta * zetal;
}

/**
 * Calculates 2-loop bosonic non-Yukawa contributions.
 *
 * Eq (71), arxiv:1607:06292
 */
double amu2L_B_nonYuk()
{
   // @todo(alex) implementation missing
   return 0.0;
}

/**
 * Calculates 2-loop bosonic Yukawa contributions.
 *
 * Eq (52), arxiv:1607:06292
 */
double amu2L_B_Yuk()
{
   // @todo(alex) pass parameters to function
   const auto tb = 1.0;
   const auto zetal = 0.0;
   const auto lambda5 = 0.0;
   const auto eta = 0.0;
   const auto al = 1.0/137;
   const auto mhSM2 = 125.0;
   const auto mH2 = 0.0;
   const auto mHp2 = 0.0;
   const auto mw2 = 82.0;
   const auto mz2 = 91.0;
   const auto mm2 = 0.1;

   const auto cw2 = mw2/mz2;
   const auto xhSM = mhSM2/mz2;
   const auto xHp = mHp2/mz2;
   const auto xH = mH2/mz2;

   // Eq.(91), arxiv:1607.06292
   const auto a000 = fb(xhSM, xHp, al, cw2)*Fm0(xhSM, xHp, al, cw2, mm2, mz2);
   // Eq.(92), arxiv:1607.06292
   const auto a0z0 = -fb(xH, 0.0, al, cw2)*(
      Fm0(xH, xHp, al, cw2, mm2, mz2) + Fmp(xH, xHp, al, cw2, mm2, mz2));
   // Eq.(93), arxiv:1607.06292
   const auto a500 = Fm0(xhSM, xHp, al, cw2, mm2, mz2);
   // Eq.(94), arxiv:1607.06292
   const auto a5z0 = -0.5*(
      Fm0(xH, xHp, al, cw2, mm2, mz2) + Fmp(xH, xHp, al, cw2, mm2, mz2));
   // Eq.(95), arxiv:1607.06292
   const auto a001 =
      + fb(xH, 0.0, al, cw2)*Fm0(xH, xHp, al, cw2, mm2, mz2)
      - fb(xhSM, 0.0, al, cw2)*Fm0(xhSM, xHp, al, cw2, mm2, mz2);
   // Eq.(96), arxiv:1607.06292
   const auto a0z1 = 0.0; // @todo(alex) implementation missing
   // Eq.(97), arxiv:1607.06292
   const auto a501 =
      + Fm0(xH, xHp, al, cw2, mm2, mz2)/2
      - Fm0(xhSM, xHp, al, cw2, mm2, mz2)/2;
   // Eq.(98), arxiv:1607.06292
   const auto a5z1 =
      - Fm0(xH, xHp, al, cw2, mm2, mz2)
      - Fmp(xH, xHp, al, cw2, mm2, mz2)
      + Fm0(xhSM, xHp, al, cw2, mm2, mz2)
      + Fmp(xhSM, xHp, al, cw2, mm2, mz2);

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

   return res;
}

} // namespace general_thdm

} // namespace gm2calc
