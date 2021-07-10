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

#include "gm2calc/SM.hpp"
#include "gm2_numerics.hpp"
#include "gm2_constants.h"
#include <cmath>

namespace gm2calc {

namespace {

const double pi = 3.1415926535897932;

} // anonymous namespace

SM::SM()
   : alpha_em_0(GM2CALC_ALPHA_EM_THOMPSON)
   , alpha_em_mz(GM2CALC_ALPHA_EM_MZ)
   , mh(GM2CALC_MH)
   , mw(GM2CALC_MW)
   , mz(GM2CALC_MZ)
   , mu((Eigen::Matrix<double,3,1>() << GM2CALC_MU, GM2CALC_MC, GM2CALC_MT).finished())
   , md((Eigen::Matrix<double,3,1>() << GM2CALC_MD, GM2CALC_MS, GM2CALC_MBMB).finished())
   , ml((Eigen::Matrix<double,3,1>() << GM2CALC_ME, GM2CALC_MM, GM2CALC_ML).finished())
{
}

double SM::get_e_0() const
{
   return std::sqrt(alpha_em_0*4*pi);
}

double SM::get_e_mz() const
{
   return std::sqrt(alpha_em_mz*4*pi);
}

double SM::get_gY() const
{
   return get_e_mz()/get_cw();
}

double SM::get_g2() const
{
   return get_e_mz()/get_sw();
}

double SM::get_cw() const
{
   return std::abs(mw/mz);
}

double SM::get_sw() const
{
   return std::sqrt(1 - sqr(get_cw()));
}

double SM::get_v() const
{
   return 2*mw/get_g2();
}


} // namespace gm2calc
