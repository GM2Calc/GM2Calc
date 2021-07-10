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
#include <cmath>

namespace gm2calc {

namespace {

const double pi = 3.1415926535897932;

} // anonymous namespace

double SM::get_e() const
{
   return std::sqrt(alpha_em*4*pi);
}

double SM::get_gY() const
{
   return get_e()/get_cw();
}

double SM::get_g2() const
{
   return get_e()/get_sw();
}

double SM::get_cw() const
{
   return std::abs(mw/mz);
}

double SM::get_sw() const
{
   return std::sqrt(1 - sqr(get_cw()));
}


} // namespace gm2calc
