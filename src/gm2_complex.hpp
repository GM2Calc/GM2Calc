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

#ifndef GM2_COMPLEX_HPP
#define GM2_COMPLEX_HPP

#include <cmath>
#include <complex>

namespace gm2calc {

template <typename T>
struct Complex {
   constexpr Complex(T re_ = T{}, T im_ = T{}) : re(re_), im(im_) {}
   operator std::complex<T>() const noexcept { return std::complex<T>(re, im); }
   T re{};
   T im{};
};

template <typename T>
constexpr T arg(const Complex<T>& z) noexcept
{
   return std::atan2(z.im, z.re);
}

template <typename T>
constexpr Complex<T> conj(const Complex<T>& z) noexcept
{
   return { z.re, -z.im };
}

template <typename T>
Complex<T> log(const Complex<T>& z) noexcept
{
   T a = arg(z);

   if (z.im == T(0) && a < T(0)) {
      a = -a;
   }

   return { 0.5*std::log(norm_sqr(z)), a };
}

template <typename T>
constexpr T norm_sqr(const Complex<T>& z) noexcept
{
   return z.re*z.re + z.im*z.im;
}

template <typename T>
constexpr Complex<T> operator+(const Complex<T>& a, const Complex<T>& b) noexcept
{
   return { a.re + b.re, a.im + b.im };
}

template <typename T>
constexpr Complex<T> operator+(const Complex<T>& z, T x) noexcept
{
   return { z.re + x, z.im };
}

template <typename T>
constexpr Complex<T> operator+(T x, const Complex<T>& z) noexcept
{
   return { x + z.re, z.im };
}

template <typename T>
constexpr Complex<T> operator-(const Complex<T>& a, const Complex<T>& b) noexcept
{
   return { a.re - b.re, a.im - b.im };
}

template <typename T>
constexpr Complex<T> operator-(T x, const Complex<T>& z) noexcept
{
   return { x - z.re, -z.im };
}

template <typename T>
constexpr Complex<T> operator-(const Complex<T>& z, T x) noexcept
{
   return { z.re - x, z.im };
}

template <typename T>
constexpr Complex<T> operator-(const Complex<T>& z) noexcept
{
   return { -z.re, -z.im };
}

template <typename T>
constexpr Complex<T> operator*(const Complex<T>& a, const Complex<T>& b) noexcept
{
   return { a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re };
}

template <typename T>
constexpr Complex<T> operator*(T x, const Complex<T>& z) noexcept
{
   return { x*z.re, x*z.im };
}

template <typename T>
constexpr Complex<T> operator*(const Complex<T>& z, T x) noexcept
{
   return x*z;
}

template <typename T>
constexpr Complex<T> operator/(T x, const Complex<T>& z) noexcept
{
   return x*conj(z)/norm_sqr(z);
}

template <typename T>
constexpr Complex<T> operator/(const Complex<T>& z, T x) noexcept
{
   return { z.re/x, z.im/x };
}

} // namespace gm2calc

#endif
