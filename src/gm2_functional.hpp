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

#ifndef GM2_FUNCTIONAL_HPP
#define GM2_FUNCTIONAL_HPP

#include <cmath>
#include <Eigen/Core>

namespace gm2calc {
namespace functional {

template<class Real, int N>
struct Abs_less {
    explicit Abs_less(const Eigen::Array<Real, N, 1>& w_) : w(w_) {}
    bool operator() (int i, int j) { return std::abs(w[i]) < std::abs(w[j]); }
    const Eigen::Array<Real, N, 1>& w;
};

template <class T>
struct Is_not_finite {
   bool operator()(T x) { return !std::isfinite(x); }
};

} // namespace functional
} // namespace gm2calc

#endif
