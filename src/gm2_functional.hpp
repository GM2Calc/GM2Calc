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

#ifndef GM2_COMPARE_H
#define GM2_COMPARE_H

#include <cmath>
#include <Eigen/Core>

namespace gm2calc {
namespace functional {

template<class Real, int N>
struct Abs_less {
    Abs_less(const Eigen::Array<Real, N, 1>& w_) : w(w_) {}
    bool operator() (int i, int j) { return std::abs(w[i]) < std::abs(w[j]); }
    const Eigen::Array<Real, N, 1>& w;
};

template<class Real>
struct Flip_sign {
    std::complex<Real> operator() (const std::complex<Real>& z) const {
	return z.real() < 0 ? std::complex<Real>(0,1) :
	    std::complex<Real>(1,0);
    }
};

template <class T>
struct Is_not_finite {
   bool operator()(T x) { return !std::isfinite(x); }
};

template<class Real, int N>
struct Less {
    Less(const Eigen::Array<Real, N, 1>& s_) : s(s_) {}
    bool operator() (int i, int j) { return s[i] < s[j]; }
    const Eigen::Array<Real, N, 1>& s;
};

template<class Real>
struct Rephase {
    std::complex<Real> operator() (const std::complex<Real>& z) const
	{ return std::polar(Real(1), std::arg(z)/2); }
};

} // namespace functional
} // namespace gm2calc

#endif
