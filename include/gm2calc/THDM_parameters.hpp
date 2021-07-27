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

#ifndef THDM_PARAMETERS_H
#define THDM_PARAMETERS_H

#include <complex>
#include <iosfwd>
#include <Eigen/Core>

namespace gm2calc {

/**
 * @class THDM_parameters
 * @brief Contains the parameters of the THDM model
 */
class THDM_parameters {
public:
   void print(std::ostream&) const;

   void set_g1(double g1_) { g1 = g1_; }
   void set_g2(double g2_) { g2 = g2_; }
   void set_g3(double g3_) { g3 = g3_; }
   void set_lambda1(double lambda1_) { lambda1 = lambda1_; }
   void set_lambda2(double lambda2_) { lambda2 = lambda2_; }
   void set_lambda3(double lambda3_) { lambda3 = lambda3_; }
   void set_lambda4(double lambda4_) { lambda4 = lambda4_; }
   void set_lambda5(double lambda5_) { lambda5 = lambda5_; }
   void set_lambda6(double lambda6_) { lambda6 = lambda6_; }
   void set_lambda7(double lambda7_) { lambda7 = lambda7_; }
   void set_m122(double m122_) { m122 = m122_; }
   void set_m112(double m112_) { m112 = m112_; }
   void set_m222(double m222_) { m222 = m222_; }
   void set_v1(double v1_) { v1 = v1_; }
   void set_v2(double v2_) { v2 = v2_; }
   void set_Yu(const Eigen::Matrix<std::complex<double>,3,3>& Yu_) { Yu = Yu_; }
   void set_Yu(int i, int k, const std::complex<double>& value) { Yu(i,k) = value; }
   void set_Xu(const Eigen::Matrix<std::complex<double>,3,3>& Xu_) { Xu = Xu_; }
   void set_Xu(int i, int k, const std::complex<double>& value) { Xu(i,k) = value; }
   void set_Yd(const Eigen::Matrix<std::complex<double>,3,3>& Yd_) { Yd = Yd_; }
   void set_Yd(int i, int k, const std::complex<double>& value) { Yd(i,k) = value; }
   void set_Yl(const Eigen::Matrix<std::complex<double>,3,3>& Yl_) { Yl = Yl_; }
   void set_Yl(int i, int k, const std::complex<double>& value) { Yl(i,k) = value; }
   void set_Xd(const Eigen::Matrix<std::complex<double>,3,3>& Xd_) { Xd = Xd_; }
   void set_Xd(int i, int k, const std::complex<double>& value) { Xd(i,k) = value; }
   void set_Xl(const Eigen::Matrix<std::complex<double>,3,3>& Xl_) { Xl = Xl_; }
   void set_Xl(int i, int k, const std::complex<double>& value) { Xl(i,k) = value; }

   double get_m122() const { return m122; }
   double get_m112() const { return m112; }
   double get_m222() const { return m222; }
   double get_v1() const { return v1; }
   double get_v2() const { return v2; }
   double get_g1() const { return g1; }
   double get_g2() const { return g2; }
   double get_g3() const { return g3; }
   double get_lambda1() const { return lambda1; }
   double get_lambda2() const { return lambda2; }
   double get_lambda3() const { return lambda3; }
   double get_lambda4() const { return lambda4; }
   double get_lambda5() const { return lambda5; }
   double get_lambda6() const { return lambda6; }
   double get_lambda7() const { return lambda7; }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Yu() const { return Yu; }
   std::complex<double> get_Yu(int i, int k) const { return Yu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Xu() const { return Xu; }
   std::complex<double> get_Xu(int i, int k) const { return Xu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Yd() const { return Yd; }
   std::complex<double> get_Yd(int i, int k) const { return Yd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Yl() const { return Yl; }
   std::complex<double> get_Yl(int i, int k) const { return Yl(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Xd() const { return Xd; }
   std::complex<double> get_Xd(int i, int k) const { return Xd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Xl() const { return Xl; }
   std::complex<double> get_Xl(int i, int k) const { return Xl(i,k); }

protected:
   double g1{0.0};
   double g2{0.0};
   double g3{0.0};
   double lambda6{0.0};
   double lambda5{0.0};
   double lambda7{0.0};
   double lambda1{0.0};
   double lambda4{0.0};
   double lambda3{0.0};
   double lambda2{0.0};
   double m122{0.0};
   double m112{0.0};
   double m222{0.0};
   double v1{0.0};
   double v2{0.0};
   Eigen::Matrix<std::complex<double>,3,3> Yu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Xu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Yd{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Yl{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Xd{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Xl{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
};

std::ostream& operator<<(std::ostream&, const THDM_parameters&);

} // namespace gm2calc

#endif