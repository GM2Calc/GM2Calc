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
 * \file gm2_2loop_F.cpp
 *
 * Contains functions necessary to calculate the fermionic THDM
 * contributions for g-2 at the 2-loop level.
 */

namespace gm2calc {

namespace {

struct F_parameters {
   double alpha{}; ///< alpha_em
   double mm{};    ///< muon mass
   double mw{};    ///< W boson mass
   double mz{};    ///< Z boson mass
   double qf{};    ///< electromagnetic charge of fermion f
   double ql{};    ///< electromagnetic charge of fermion l
   double t3f{};   ///< SU(2)_L charge of ferimon f
   double t3l{};   ///< SU(2)_L charge of ferimon l
   double nc{};    ///< number of colors of fermion f
};

const double pi = 3.1415926535897932;

double sqr(double x) noexcept { return x*x; }

/// Eq (56), arxiv:1607.06292, S = h or H
double FS(double ms2, double mf2)
{
   return -2.0 + std::log(ms2/mf2)
      - (ms2 - 2*mf2)/ms2 * Phi(ms2, mf2, mf2)/(ms2 - 4*mf2);
}

/// Eq (57), arxiv:1607.06292, S = A
double FA(double ms2, double mf2)
{
   return Phi(ms2, mf2, mf2)/(ms2 - 4*mf2);
}

/// Eq (54), arxiv:1607.06292, S = h or H
double fSgamma(double ms2, double mf2, const F_parameters& pars) noexcept
{
   const double al2 = sqr(pars.alpha);
   const double mm2 = sqr(pars.mm);
   const double mw2 = sqr(pars.mw);
   const double mz2 = sqr(pars.mz);
   const double sw2 = 1.0 - mw2/mz2;
   const double qf2 = sqr(pars.qf);
   const double nc = pars.nc;

   return al2*mm2/(4*sqr(pi)*mw2*sw2) * qf2*nc * mf2/ms2 * FS(ms2, mf2);
}

/// Eq (55), arxiv:1607.06292, S = h or H
double fSZ(double ms2, double mf2, const F_parameters& pars) noexcept
{
   const double al2 = sqr(pars.alpha);
   const double mm2 = sqr(pars.mm);
   const double mw2 = sqr(pars.mw);
   const double mz2 = sqr(pars.mz);
   const double cw2 = mw2/mz2;
   const double sw2 = 1.0 - cw2;
   const double qf = pars.qf;
   const double ql = pars.ql;
   const double nc = pars.nc;
   const double gvf = 0.5*pars.t3f - qf*sw2;
   const double gvl = 0.5*pars.t3l - ql*sw2;

   return al2*mm2/(4*sqr(pi)*mw2*sw2) * (-nc*qf*gvl*gvf)/(sw2*cw2)
      * mf2/(ms2 - mz2) * (FS(ms2, mf2) - FS(mz2, mf2));
}

/// Eq (53), arxiv:1607.06292, S = h or H
double ffS(double x, double y, const F_parameters& pars) noexcept
{
   return fSgamma(x, y, pars) + fSZ(x, y, pars);
}

/// Eq (60), arxiv:1607.06292
double FlHp(double ms2, double mf2) noexcept
{
   const double xl = mf2/ms2;

   return xl + xl*(xl - 1.0)*(dilog(1.0 - 1.0/xl) - sqr(pi)/6)
      + (xl - 0.5)*std::log(xl);
}

/// Eq (61), arxiv:1607.06292
double FdHp(double ms2, double md2, double mu2, double qd, double qu) noexcept
{
   const double xu = mu2/ms2;
   const double xd = md2/ms2;
   const double sqrt_xu = std::sqrt(xu);
   const double sqrt_xd = std::sqrt(xd);
   const double y = sqr(xu - xd) - 2*(xu + xd) + 1.0;
   const double s = 0.25*(qu + qd);
   const double c = sqr(xu - xd) - qu*xu + qd*xd;
   const double cbar = (xu - qu)*xu - (xd + qd)*xd;

   return -(xu - xd) + (cbar/y - c*(xu - xd)/y) * Phi(sqrt_xd, sqrt_xu, 1.0)
      + c*(dilog(1.0 - xd/xu) - 0.5*std::log(xu)*std::log(xd/xu)*Phi(sqrt_xd, sqrt_xu, 1.0))
      + (s + xd)*std::log(xd) + (s - xu)*std::log(xu);
}

/// Eq (62), arxiv:1607.06292
double FuHp(double ms2, double md2, double mu2, double qd, double qu) noexcept
{
   const double xu = mu2/ms2;
   const double xd = md2/ms2;
   const double sqrt_xu = std::sqrt(xu);
   const double sqrt_xd = std::sqrt(xd);
   const double y = sqr(xu - xd) - 2*(xu + xd) + 1.0;

   return FdHp(ms2, md2, mu2, 2.0 + qd, 2.0 + qu)
      - 4.0/3.0*(xu - xd - 1.0)/y*Phi(sqrt_xd, sqrt_xu, 1.0)
      - 1.0/3.0*(sqr(std::log(xd)) - sqr(std::log(xu)));
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = l
double flHp(double ms2, double mf2, const F_parameters& pars) noexcept
{
   const double al2 = sqr(pars.alpha);
   const double mm2 = sqr(pars.mm);
   const double mw2 = sqr(pars.mw);
   const double mz2 = sqr(pars.mz);
   const double cw2 = mw2/mz2;
   const double sw2 = 1.0 - cw2;
   const double sw4 = sqr(sw2);
   const double nc = pars.nc;

   return al2*mm2/(32*sqr(pi)*mw2*sw4) * nc*mf2/(ms2 - mw2)
      * (FlHp(ms2, mf2) - FlHp(mw2, mf2));
}

} // anonymous namespace

/**
 * \fn amu2L_F
 *
 * Calculates 2-loop fermionic contributions.
 *
 * Eq (63), arxiv:1607:06292
 */
double amu2L_F()
{
   // @todo(alex) implementation missing
   const double res = 0.0;
   return res;
}

} // namespace gm2calc
