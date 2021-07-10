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

#ifndef GM2_SM_HPP
#define GM2_SM_HPP

#include <Eigen/Core>

namespace gm2calc {

class SM {
public:
   void set_alpha_em(double a) noexcept { alpha_em = a; }
   void set_mw(double m) noexcept { mw = m; }
   void set_mz(double m) noexcept { mz = m; }
   void set_mh(double m) noexcept { mh = m; }
   void set_mu(const Eigen::Matrix<double,3,1>& m) noexcept { mu = m; }
   void set_md(const Eigen::Matrix<double,3,1>& m) noexcept { md = m; }
   void set_ml(const Eigen::Matrix<double,3,1>& m) noexcept { ml = m; }
   void set_mu(int i, double m) noexcept { mu(i) = m; }
   void set_md(int i, double m) noexcept { md(i) = m; }
   void set_ml(int i, double m) noexcept { ml(i) = m; }

   double get_alpha_em() const { return alpha_em; }
   double get_mw() const { return mw; }
   double get_mz() const { return mz; }
   double get_mh() const { return mh; }
   const Eigen::Matrix<double,3,1>& get_mu() const { return mu; }
   const Eigen::Matrix<double,3,1>& get_md() const { return md; }
   const Eigen::Matrix<double,3,1>& get_ml() const { return ml; }
   double get_mu(int i) const { return mu(i); }
   double get_md(int i) const { return md(i); }
   double get_ml(int i) const { return ml(i); }
   double get_e() const;
   double get_gY() const;
   double get_g2() const;
   double get_cw() const;
   double get_sw() const;

private:
   double alpha_em{0.0};
   double mw{0.0};
   double mz{0.0};
   double mh{0.0};
   Eigen::Matrix<double,3,1> mu{Eigen::Matrix<double,3,1>::Zero()}; ///< up-type quark masses
   Eigen::Matrix<double,3,1> md{Eigen::Matrix<double,3,1>::Zero()}; ///< down-type quark masses
   Eigen::Matrix<double,3,1> ml{Eigen::Matrix<double,3,1>::Zero()}; ///< down-type lepton masses
};

} // namespace gm2calc

#endif
