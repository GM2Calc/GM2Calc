// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef GM2_NUMERICS_HPP
#define GM2_NUMERICS_HPP

#include <cmath>
#include <limits>
#include <cstddef>

namespace gm2calc {

/// returns number squared
template <typename T> T sqr(T x) noexcept { return x*x; }

/// returns square root of absolute of number
double abs_sqrt(double);

/// returns sign of real number
int sign(double) noexcept;

/// returns square root of absolute of number, times sign
double signed_abs_sqrt(double);

/// returns square of number, times sign
double signed_sqr(double);

template <typename T>
bool is_zero(T a, T eps)
{
   return std::fabs(a) < eps;
}

template <typename T>
bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
{
   return is_zero(a - b, prec);
}

template <typename T>
bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon())
{
   if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
      return true;

   if (std::fabs(a) < std::numeric_limits<T>::epsilon() ||
       std::fabs(b) < std::numeric_limits<T>::epsilon())
      return false;

   return std::fabs((a - b)/a) < prec;
}

} // namespace gm2calc

#endif
