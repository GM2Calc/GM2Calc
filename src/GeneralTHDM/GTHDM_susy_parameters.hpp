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

#ifndef GTHDM_susy_parameters_H
#define GTHDM_susy_parameters_H

#include <iosfwd>

#include <Eigen/Core>

namespace gm2calc {

/**
 * @class GTHDM_susy_parameters
 * @brief contains dimensionless parameters of the GTHDM model
 *
 * dimensionless parameters are: Gauge couplings, Yukawa couplings and
 * scalar couplings.  In addition, this class stores the current
 * renormalization scale.
 */
class GTHDM_susy_parameters {
public:
   virtual ~GTHDM_susy_parameters() = default;

   virtual void print(std::ostream&) const;

   void set_scale(double s) { scale = s; }
   double get_scale() const { return scale; }

   void set_g1(double g1_) { g1 = g1_; }
   void set_g2(double g2_) { g2 = g2_; }
   void set_g3(double g3_) { g3 = g3_; }
   void set_Lambda1(double Lambda1_) { Lambda1 = Lambda1_; }
   void set_Lambda2(double Lambda2_) { Lambda2 = Lambda2_; }
   void set_Lambda3(double Lambda3_) { Lambda3 = Lambda3_; }
   void set_Lambda4(double Lambda4_) { Lambda4 = Lambda4_; }
   void set_Lambda5(double Lambda5_) { Lambda5 = Lambda5_; }
   void set_Lambda6(double Lambda6_) { Lambda6 = Lambda6_; }
   void set_Lambda7(double Lambda7_) { Lambda7 = Lambda7_; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) { Yu = Yu_; }
   void set_Yu(int i, int k, const double& value) { Yu(i,k) = value; }
   void set_Xu(const Eigen::Matrix<double,3,3>& Xu_) { Xu = Xu_; }
   void set_Xu(int i, int k, const double& value) { Xu(i,k) = value; }
   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) { Yd = Yd_; }
   void set_Yd(int i, int k, const double& value) { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) { Ye = Ye_; }
   void set_Ye(int i, int k, const double& value) { Ye(i,k) = value; }
   void set_Xd(const Eigen::Matrix<double,3,3>& Xd_) { Xd = Xd_; }
   void set_Xd(int i, int k, const double& value) { Xd(i,k) = value; }
   void set_Xe(const Eigen::Matrix<double,3,3>& Xe_) { Xe = Xe_; }
   void set_Xe(int i, int k, const double& value) { Xe(i,k) = value; }

   double get_g1() const { return g1; }
   double get_g2() const { return g2; }
   double get_g3() const { return g3; }
   double get_Lambda1() const { return Lambda1; }
   double get_Lambda2() const { return Lambda2; }
   double get_Lambda3() const { return Lambda3; }
   double get_Lambda4() const { return Lambda4; }
   double get_Lambda5() const { return Lambda5; }
   double get_Lambda6() const { return Lambda6; }
   double get_Lambda7() const { return Lambda7; }
   const Eigen::Matrix<double,3,3>& get_Yu() const { return Yu; }
   double get_Yu(int i, int k) const { return Yu(i,k); }
   const Eigen::Matrix<double,3,3>& get_Xu() const { return Xu; }
   double get_Xu(int i, int k) const { return Xu(i,k); }
   const Eigen::Matrix<double,3,3>& get_Yd() const { return Yd; }
   double get_Yd(int i, int k) const { return Yd(i,k); }
   const Eigen::Matrix<double,3,3>& get_Ye() const { return Ye; }
   double get_Ye(int i, int k) const { return Ye(i,k); }
   const Eigen::Matrix<double,3,3>& get_Xd() const { return Xd; }
   double get_Xd(int i, int k) const { return Xd(i,k); }
   const Eigen::Matrix<double,3,3>& get_Xe() const { return Xe; }
   double get_Xe(int i, int k) const { return Xe(i,k); }

protected:
   double scale{0.0};
   double g1{0.0};
   double g2{0.0};
   double g3{0.0};
   double Lambda6{0.0};
   double Lambda5{0.0};
   double Lambda7{0.0};
   double Lambda1{0.0};
   double Lambda4{0.0};
   double Lambda3{0.0};
   double Lambda2{0.0};
   Eigen::Matrix<double,3,3> Yu{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Xu{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Yd{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Ye{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Xd{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Xe{Eigen::Matrix<double,3,3>::Zero()};
};

std::ostream& operator<<(std::ostream&, const GTHDM_susy_parameters&);

} // namespace gm2calc

#endif
