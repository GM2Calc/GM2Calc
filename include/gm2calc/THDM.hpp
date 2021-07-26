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

#ifndef GM2_THDM_HPP
#define GM2_THDM_HPP

#include "THDM_mass_eigenstates.hpp"
#include "SM.hpp"

#include <iosfwd>

#include <Eigen/Core>

namespace gm2calc {

/**
 * @class THDM
 * @brief Contains routines to determine the THDM parameters
 */
class THDM : private THDM_mass_eigenstates {
public:
   enum class Yukawa_scheme {
      type_1, type_2, type_X, type_Y, general
   };

   struct General_basis {
      Yukawa_scheme yukawa_scheme{Yukawa_scheme::type_2};
      double lambda1{0.0};
      double lambda2{0.0};
      double lambda3{0.0};
      double lambda4{0.0};
      double lambda5{0.0};
      double lambda6{0.0};
      double lambda7{0.0};
      double tan_beta{0.0};
      double m122{0.0};
      Eigen::Matrix<double,3,3> Xu{Eigen::Matrix<double,3,3>::Zero()};
      Eigen::Matrix<double,3,3> Xd{Eigen::Matrix<double,3,3>::Zero()};
      Eigen::Matrix<double,3,3> Xl{Eigen::Matrix<double,3,3>::Zero()};
   };

   struct Physical_basis {
      Yukawa_scheme yukawa_scheme{Yukawa_scheme::type_2};
      double mh{0.0};
      double mH{0.0};
      double mA{0.0};
      double mHp{0.0};
      double sin_beta_minus_alpha{0.0};
      double lambda6{0.0};
      double lambda7{0.0};
      double tan_beta{0.0};
      double m122{0.0};
      Eigen::Matrix<double,3,3> Xu{Eigen::Matrix<double,3,3>::Zero()};
      Eigen::Matrix<double,3,3> Xd{Eigen::Matrix<double,3,3>::Zero()};
      Eigen::Matrix<double,3,3> Xl{Eigen::Matrix<double,3,3>::Zero()};
   };

   THDM(const General_basis&, const SM& sm_ = SM{});
   THDM(const Physical_basis&, const SM& sm_ = SM{});

   double get_zeta_u() const;
   double get_zeta_d() const;
   double get_zeta_l() const;

   Eigen::Matrix<double,3,3> get_xi_u() const;
   Eigen::Matrix<double,3,3> get_xi_d() const;
   Eigen::Matrix<double,3,3> get_xi_l() const;

   Eigen::Matrix<double,3,3> get_yuh() const;
   Eigen::Matrix<double,3,3> get_yuH() const;
   Eigen::Matrix<double,3,3> get_yuA() const;
   Eigen::Matrix<std::complex<double>,3,3> get_yuHp() const;

   Eigen::Matrix<double,3,3> get_ydh() const;
   Eigen::Matrix<double,3,3> get_ydH() const;
   Eigen::Matrix<double,3,3> get_ydA() const;
   Eigen::Matrix<std::complex<double>,3,3> get_ydHp() const;

   Eigen::Matrix<double,3,3> get_ylh() const;
   Eigen::Matrix<double,3,3> get_ylH() const;
   Eigen::Matrix<double,3,3> get_ylA() const;
   Eigen::Matrix<double,3,3> get_ylHp() const;

   Yukawa_scheme get_yukawa_scheme() const { return yukawa_scheme; }

   const SM& get_sm() const { return sm; }

   void set_tan_beta(double);

   using THDM_mass_eigenstates::get_alpha_em;
   using THDM_mass_eigenstates::get_alpha_h;
   using THDM_mass_eigenstates::get_beta;
   using THDM_mass_eigenstates::get_eta;
   using THDM_mass_eigenstates::get_tan_beta;
   using THDM_mass_eigenstates::get_v;
   using THDM_mass_eigenstates::get_v_sqr;
   using THDM_mass_eigenstates::get_lambda1;
   using THDM_mass_eigenstates::get_lambda2;
   using THDM_mass_eigenstates::get_lambda3;
   using THDM_mass_eigenstates::get_lambda4;
   using THDM_mass_eigenstates::get_lambda5;
   using THDM_mass_eigenstates::get_lambda6;
   using THDM_mass_eigenstates::get_lambda7;
   using THDM_mass_eigenstates::get_LambdaFive;
   using THDM_mass_eigenstates::get_m122;
   using THDM_mass_eigenstates::get_g1;
   using THDM_mass_eigenstates::get_g2;
   using THDM_mass_eigenstates::get_Yu;
   using THDM_mass_eigenstates::get_Yd;
   using THDM_mass_eigenstates::get_Yl;
   using THDM_mass_eigenstates::get_Xu;
   using THDM_mass_eigenstates::get_Xd;
   using THDM_mass_eigenstates::get_Xl;
   using THDM_mass_eigenstates::get_v1;
   using THDM_mass_eigenstates::get_v2;
   using THDM_mass_eigenstates::get_Mhh;
   using THDM_mass_eigenstates::get_MAh;
   using THDM_mass_eigenstates::get_MHm;
   using THDM_mass_eigenstates::get_MFu;
   using THDM_mass_eigenstates::get_MFd;
   using THDM_mass_eigenstates::get_MFv;
   using THDM_mass_eigenstates::get_MFe;
   using THDM_mass_eigenstates::get_MVG;
   using THDM_mass_eigenstates::get_MVP;
   using THDM_mass_eigenstates::get_MVWm;
   using THDM_mass_eigenstates::get_MVZ;

   friend std::ostream& operator<<(std::ostream&, const THDM&);

private:
   SM sm{};
   Yukawa_scheme yukawa_scheme{Yukawa_scheme::type_2};

   void init_gauge_couplings();
   void init_yukawas();
   void set_basis(const General_basis&);
   void set_basis(const Physical_basis&);
};

/// streaming operator
std::ostream& operator<<(std::ostream&, const THDM&);

} // namespace gm2calc

#endif
