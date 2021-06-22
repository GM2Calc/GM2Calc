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

#include <cmath>

/**
 * \file gm2_2loop_F.cpp
 *
 * Contains functions necessary to calculate the fermionic THDM
 * contributions for g-2 at the 2-loop level.
 */

namespace gm2calc {

namespace general_thdm {

namespace {

struct F_sm_pars {
   double alpha2; ///< squared alpha_em
   double mm2;    ///< squared muon mass
   double mw2;    ///< squared W boson mass
   double mz2;    ///< squared Z boson mass
};

struct F_neut_pars {
   double mf2;   ///< squared mass of fermion f
   double qf;    ///< electromagnetic charge of fermion f
   double ql;    ///< electromagnetic charge of fermion l
   double t3f;   ///< SU(2)_L charge of ferimon f
   double t3l;   ///< SU(2)_L charge of ferimon l
   double nc;    ///< number of colors of fermion f
};

struct F_char_pars {
   double mf2;   ///< squared fermion mass
   double md2;   ///< squared mass of up-type fermion
   double mu2;   ///< squared mass of down-type fermion
   double qd;    ///< electromagnetic charge of up-type fermion
   double qu;    ///< electromagnetic charge of down-type fermion
   double nc;    ///< number of colors
};

const double pi = 3.1415926535897932;
const double q_u =  2.0/3.0;  ///< electric charge of up-type quarks
const double q_d = -1.0/3.0;  ///< electric charge of down-type quarks
const double q_v =  0.0;      ///< electric charge of up-type leptons
const double q_l = -1.0;      ///< electric charge of down-type leptons
const double t3_u = +1.0;     ///< SU(2)_L charge of up-type quark @todo(alex) check convention
const double t3_d = -1.0;     ///< SU(2)_L charge of down-type quark @todo(alex) check convention
const double t3_l = -1.0;     ///< SU(2)_L charge of charged lepton @todo(alex) check convention

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
template <typename F>
double fSgamma(double ms2, const F_neut_pars& pars, const F_sm_pars& sm, F FS) noexcept
{
   const double al2 = sm.alpha2;
   const double mm2 = sm.mm2;
   const double mw2 = sm.mw2;
   const double mz2 = sm.mz2;
   const double sw2 = 1.0 - mw2/mz2;
   const double mf2 = pars.mf2;
   const double qf2 = sqr(pars.qf);
   const double nc = pars.nc;

   return al2*mm2/(4*sqr(pi)*mw2*sw2) * qf2*nc * mf2/ms2 * FS(ms2, mf2);
}

/// Eq (55), arxiv:1607.06292, S = h or H
template <typename F>
double fSZ(double ms2, const F_neut_pars& pars, const F_sm_pars& sm, F FS) noexcept
{
   const double al2 = sm.alpha2;
   const double mm2 = sm.mm2;
   const double mw2 = sm.mw2;
   const double mz2 = sm.mz2;
   const double cw2 = mw2/mz2;
   const double sw2 = 1.0 - cw2;
   const double mf2 = pars.mf2;
   const double qf = pars.qf;
   const double ql = pars.ql;
   const double nc = pars.nc;
   const double gvf = 0.5*pars.t3f - qf*sw2;
   const double gvl = 0.5*pars.t3l - ql*sw2;

   return al2*mm2/(4*sqr(pi)*mw2*sw2) * (-nc*qf*gvl*gvf)/(sw2*cw2)
      * mf2/(ms2 - mz2) * (FS(ms2, mf2) - FS(mz2, mf2));
}

/// Eq (53), arxiv:1607.06292, S = h or H
template <typename F>
double ffS(double ms2, const F_neut_pars& pars, const F_sm_pars& sm, F FS) noexcept
{
   return fSgamma(ms2, pars, sm, FS) + fSZ(ms2, pars, sm, FS);
}

/// Eq (60), arxiv:1607.06292
double FlHp(double ms2, double mf2) noexcept
{
   const double xl = mf2/ms2;

   return xl + xl*(xl - 1.0)*(dilog(1.0 - 1.0/xl) - sqr(pi)/6)
      + (xl - 0.5)*std::log(xl);
}

/**
 * Eq (61), arxiv:1607.06292
 *
 * @note There is a misprint in Eq (61), arxiv:1607.06292v2: There
 * should be no Phi function in the 2nd line of (61).
 *
 * @return
 */
double FdHp(double ms2, double md2, double mu2, double qd, double qu) noexcept
{
   const double xu = mu2/ms2;
   const double xd = md2/ms2;
   const double y = sqr(xu - xd) - 2*(xu + xd) + 1.0;
   const double s = 0.25*(qu + qd);
   const double c = sqr(xu - xd) - qu*xu + qd*xd;
   const double cbar = (xu - qu)*xu - (xd + qd)*xd;
   const double lxu = std::log(xu);

   return -(xu - xd) + (cbar/y - c*(xu - xd)/y) * Phi(xd, xu, 1.0)
      + c*(dilog(1.0 - xd/xu) - 0.5*lxu*std::log(xd/xu))
      + (s + xd)*std::log(xd) + (s - xu)*lxu;
}

/// Eq (62), arxiv:1607.06292
double FuHp(double ms2, double md2, double mu2, double qd, double qu) noexcept
{
   const double xu = mu2/ms2;
   const double xd = md2/ms2;
   const double y = sqr(xu - xd) - 2*(xu + xd) + 1.0;

   return FdHp(ms2, md2, mu2, 2.0 + qd, 2.0 + qu)
      - 4.0/3.0*(xu - xd - 1.0)/y*Phi(xd, xu, 1.0)
      - 1.0/3.0*(sqr(std::log(xd)) - sqr(std::log(xu)));
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = l
double flHp(double ms2, const F_char_pars& pars, const F_sm_pars& sm) noexcept
{
   const double al2 = sm.alpha2;
   const double mm2 = sm.mm2;
   const double mw2 = sm.mw2;
   const double mz2 = sm.mz2;
   const double cw2 = mw2/mz2;
   const double sw2 = 1.0 - cw2;
   const double sw4 = sqr(sw2);
   const double mf2 = pars.mf2;
   const double ml2 = pars.md2;
   const double nc = pars.nc;

   return al2*mm2/(32*sqr(pi)*mw2*sw4) * nc*mf2/(ms2 - mw2)
      * (FlHp(ms2, ml2) - FlHp(mw2, ml2));
}

/// Eq (59), arxiv:1607.06292, S = H^\pm, f = u or d
template <typename F>
double fqHp(double ms2, const F_char_pars& pars, const F_sm_pars& sm, F FfHp) noexcept
{
   const double al2 = sm.alpha2;
   const double mm2 = sm.mm2;
   const double mw2 = sm.mw2;
   const double mz2 = sm.mz2;
   const double cw2 = mw2/mz2;
   const double sw2 = 1.0 - cw2;
   const double sw4 = sqr(sw2);
   const double mf2 = pars.mf2;
   const double md2 = pars.md2;
   const double mu2 = pars.mu2;
   const double qd = pars.qd;
   const double qu = pars.qu;
   const double nc = pars.nc;

   return al2*mm2/(32*sqr(pi)*mw2*sw4) * nc*mf2/(ms2 - mw2)
      * (FfHp(ms2, md2, mu2, qd, qu) - FfHp(mw2, md2, mu2, qd, qu))
      ;
}

} // anonymous namespace

/**
 * Calculates 2-loop fermionic contributions with charged Higgs
 * bosons.
 *
 * Eq (58), arxiv:1607:06292
 */
double amu2L_F_charged(const THDM_F_parameters& thdm) noexcept
{
   const F_sm_pars sm{ sqr(thdm.alpha), sqr(thdm.mm), sqr(thdm.mw), sqr(thdm.mz) };
   const double mHp2 = sqr(thdm.mHp);

   const auto lFuHp = [] (double ms2, double md2, double mu2, double qd, double qu) {
      return FuHp(ms2, md2, mu2, qd, qu);
   };
   const auto lFdHp = [] (double ms2, double md2, double mu2, double qd, double qu) {
      return FdHp(ms2, md2, mu2, qd, qu);
   };

   double res = 0.0;

   // loop over generations
   for (int i = 0; i < 3; ++i) {
      const F_char_pars pars_u{sqr(thdm.mu(i)), sqr(thdm.md(i)), sqr(thdm.mu(i)), q_d, q_u, 3.0};
      const F_char_pars pars_d{sqr(thdm.md(i)), sqr(thdm.md(i)), sqr(thdm.mu(i)), q_d, q_u, 3.0};
      const F_char_pars pars_l{sqr(thdm.ml(i)), sqr(thdm.ml(i)), 0.0, q_l, q_v, 1.0};

      // H^\pm
      res += fqHp(mHp2, pars_u, sm, lFuHp)*thdm.yuS(i,2)*thdm.ylS(1,2);
      res += fqHp(mHp2, pars_d, sm, lFdHp)*thdm.ydS(i,2)*thdm.ylS(1,2);
      res += flHp(mHp2, pars_l, sm)*thdm.ylS(i,2)*thdm.ylS(1,2);
   }

   return res;
}

/**
 * Calculates 2-loop fermionic contributions with neutral Higgs
 * bosons.
 *
 * Eq (53), arxiv:1607:06292
 */
double amu2L_F_neutral(const THDM_F_parameters& thdm) noexcept
{
   const F_sm_pars sm{ sqr(thdm.alpha), sqr(thdm.mm), sqr(thdm.mw), sqr(thdm.mz) };
   const double mh2 = sqr(thdm.mh(0));
   const double mH2 = sqr(thdm.mh(1));
   const double mA2 = sqr(thdm.mA);
   const double mhSM2 = sqr(thdm.mhSM);

   const auto lFS = [] (double ms2, double mf2) { return FS(ms2, mf2); };
   const auto lFA = [] (double ms2, double mf2) { return FA(ms2, mf2); };

   double res = 0.0;

   // loop over generations
   for (int i = 0; i < 3; ++i) {
      const F_neut_pars pars_u{sqr(thdm.mu(i)), q_u, q_l, t3_u, t3_l, 3.0};
      const F_neut_pars pars_d{sqr(thdm.md(i)), q_d, q_l, t3_d, t3_l, 3.0};
      const F_neut_pars pars_l{sqr(thdm.ml(i)), q_l, q_l, t3_l, t3_l, 1.0};

      // h
      res += ffS(mh2, pars_u, sm, lFS)*thdm.yuS(i,0)*thdm.ylS(1,0);
      res += ffS(mh2, pars_d, sm, lFS)*thdm.ydS(i,0)*thdm.ylS(1,0);
      res += ffS(mh2, pars_l, sm, lFS)*thdm.ylS(i,0)*thdm.ylS(1,0);

      // H
      res += ffS(mH2, pars_u, sm, lFS)*thdm.yuS(i,1)*thdm.ylS(1,1);
      res += ffS(mH2, pars_d, sm, lFS)*thdm.ydS(i,1)*thdm.ylS(1,1);
      res += ffS(mH2, pars_l, sm, lFS)*thdm.ylS(i,1)*thdm.ylS(1,1);

      // A
      res += ffS(mA2, pars_u, sm, lFA)*thdm.yuS(i,2)*thdm.ylS(1,2);
      res += ffS(mA2, pars_d, sm, lFA)*thdm.ydS(i,2)*thdm.ylS(1,2);
      res += ffS(mA2, pars_l, sm, lFA)*thdm.ylS(i,2)*thdm.ylS(1,2);

      // subtract hSM
      res -= ffS(mhSM2, pars_u, sm, lFS);
      res -= ffS(mhSM2, pars_d, sm, lFS);
      res -= ffS(mhSM2, pars_l, sm, lFS);
   }

   return res;
}

/**
 * Calculates the sum of the 2-loop fermionic contributions with
 * neutral and charged Higgs bosons.
 *
 * Eq (63), arxiv:1607:06292
 */
double amu2L_F(const THDM_F_parameters& thdm) noexcept
{
   return amu2L_F_neutral(thdm) + amu2L_F_charged(thdm);
}

} // namespace general_thdm

} // namespace gm2calc
