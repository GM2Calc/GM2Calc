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

#ifndef GM2_GeneralTHDM_HPP
#define GM2_GeneralTHDM_HPP

#include "GeneralTHDM_mass_eigenstates.hpp"
#include "SM.hpp"

#include <iosfwd>

#include <Eigen/Core>

namespace gm2calc {

/**
 * @class GeneralTHDM
 * @brief Contains routines to determine the THDM parameters
 */
class GeneralTHDM : private GeneralTHDM_mass_eigenstates {
public:
   enum class Yukawa_scheme {
      type_1, type_2, type_X, type_Y
   };

   struct General_basis {
      double lambda1{0.0};
      double lambda2{0.0};
      double lambda3{0.0};
      double lambda4{0.0};
      double lambda5{0.0};
      double lambda6{0.0};
      double lambda7{0.0};
      double tan_beta{0.0};
      double M122{0.0};
      Eigen::Matrix<double,3,3> Xu{Eigen::Matrix<double,3,3>::Zero()};
      Eigen::Matrix<double,3,3> Xd{Eigen::Matrix<double,3,3>::Zero()};
      Eigen::Matrix<double,3,3> Xe{Eigen::Matrix<double,3,3>::Zero()};
   };

   struct Physical_basis {
      double mh{0.0};
      double mH{0.0};
      double mA{0.0};
      double mHp{0.0};
      double sin_beta_minus_alpha{0.0};
      double lambda6{0.0};
      double lambda7{0.0};
      double tan_beta{0.0};
      double M122{0.0};
      Eigen::Matrix<double,3,3> Xu{Eigen::Matrix<double,3,3>::Zero()};
      Eigen::Matrix<double,3,3> Xd{Eigen::Matrix<double,3,3>::Zero()};
      Eigen::Matrix<double,3,3> Xe{Eigen::Matrix<double,3,3>::Zero()};
   };

   GeneralTHDM(const General_basis&, const SM& sm_ = SM{});
   GeneralTHDM(const Physical_basis&, const SM& sm_ = SM{});

   double get_zeta_u() const;
   double get_zeta_d() const;
   double get_zeta_l() const;
   Eigen::Matrix<double,3,3> get_zeta_l_matrix() const;

   Eigen::Matrix<double,3,3> get_ylh() const;
   Eigen::Matrix<double,3,3> get_ylH() const;
   Eigen::Matrix<double,3,3> get_ylA() const;
   Eigen::Matrix<double,3,3> get_ylHp() const;

   Yukawa_scheme get_yukawa_scheme() const { return yukawa_scheme; }
   void set_yukawa_scheme(Yukawa_scheme ys) { yukawa_scheme = ys; }

   const SM& get_sm() const { return sm; }

   void set_tan_beta(double);

   using GeneralTHDM_mass_eigenstates::get_alpha_em;
   using GeneralTHDM_mass_eigenstates::get_alpha_h;
   using GeneralTHDM_mass_eigenstates::get_beta;
   using GeneralTHDM_mass_eigenstates::get_eta;
   using GeneralTHDM_mass_eigenstates::get_tan_beta;
   using GeneralTHDM_mass_eigenstates::get_v;
   using GeneralTHDM_mass_eigenstates::get_v_sqr;
   using GeneralTHDM_mass_eigenstates::get_Lambda1;
   using GeneralTHDM_mass_eigenstates::get_Lambda2;
   using GeneralTHDM_mass_eigenstates::get_Lambda3;
   using GeneralTHDM_mass_eigenstates::get_Lambda4;
   using GeneralTHDM_mass_eigenstates::get_Lambda5;
   using GeneralTHDM_mass_eigenstates::get_Lambda6;
   using GeneralTHDM_mass_eigenstates::get_Lambda7;
   using GeneralTHDM_mass_eigenstates::get_LambdaFive;
   using GeneralTHDM_mass_eigenstates::get_M122;
   using GeneralTHDM_mass_eigenstates::get_g1;
   using GeneralTHDM_mass_eigenstates::get_g2;
   using GeneralTHDM_mass_eigenstates::get_Yu;
   using GeneralTHDM_mass_eigenstates::get_Yd;
   using GeneralTHDM_mass_eigenstates::get_Ye;
   using GeneralTHDM_mass_eigenstates::get_Xu;
   using GeneralTHDM_mass_eigenstates::get_Xd;
   using GeneralTHDM_mass_eigenstates::get_Xe;
   using GeneralTHDM_mass_eigenstates::get_v1;
   using GeneralTHDM_mass_eigenstates::get_v2;
   using GeneralTHDM_mass_eigenstates::get_Mhh;
   using GeneralTHDM_mass_eigenstates::get_MAh;
   using GeneralTHDM_mass_eigenstates::get_MHm;
   using GeneralTHDM_mass_eigenstates::get_MFu;
   using GeneralTHDM_mass_eigenstates::get_MFd;
   using GeneralTHDM_mass_eigenstates::get_MFv;
   using GeneralTHDM_mass_eigenstates::get_MFe;
   using GeneralTHDM_mass_eigenstates::get_MVG;
   using GeneralTHDM_mass_eigenstates::get_MVP;
   using GeneralTHDM_mass_eigenstates::get_MVWm;
   using GeneralTHDM_mass_eigenstates::get_MVZ;

   using GeneralTHDM_mass_eigenstates::get_problems;

   friend std::ostream& operator<<(std::ostream&, const GeneralTHDM&);

private:
   Yukawa_scheme yukawa_scheme{Yukawa_scheme::type_2};
   SM sm{};

   void init_gauge_couplings();
   void init_yukawas();
   void set_basis(const General_basis&);
   void set_basis(const Physical_basis&);
};

/// streaming operator
std::ostream& operator<<(std::ostream&, const GeneralTHDM&);

} // namespace gm2calc

#endif
