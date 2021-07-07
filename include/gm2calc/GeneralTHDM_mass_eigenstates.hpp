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

/**
 * @file GeneralTHDM_mass_eigenstates.hpp
 *
 * @brief contains class for general THDM model with routines needed
 *        to solve EWSB and determine the masses and mixings.
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.3 .
 */

#ifndef GM2_GeneralTHDM_MASS_EIGENSTATES_H
#define GM2_GeneralTHDM_MASS_EIGENSTATES_H

#include "GeneralTHDM_physical.hpp"
#include "GeneralTHDM_problems.hpp"
#include "GeneralTHDM_soft_parameters.hpp"

#include <iosfwd>

#include <Eigen/Core>

namespace gm2calc {

/**
 * @class GeneralTHDM_mass_eigenstates
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
class GeneralTHDM_mass_eigenstates : public GeneralTHDM_soft_parameters
{
public:
   virtual ~GeneralTHDM_mass_eigenstates() = default;

   void print(std::ostream&) const override;

   void calculate_MSbar_masses();
   void copy_MSbar_masses_to_pole_masses();
   void do_force_output(bool);
   bool do_force_output() const;
   void reorder_MSbar_masses();
   void reorder_pole_masses();
   void set_physical(const GeneralTHDM_physical&);
   const GeneralTHDM_physical& get_physical() const;
   GeneralTHDM_physical& get_physical();
   const GeneralTHDM_problems& get_problems() const;
   GeneralTHDM_problems& get_problems();
   int solve_ewsb_tree_level();
   int solve_ewsb();

   double get_MVG() const { return MVG; }
   const Eigen::Array<double,3,1>& get_MFv() const { return MFv; }
   double get_MFv(int i) const { return MFv(i); }
   const Eigen::Array<double,2,1>& get_Mhh() const { return Mhh; }
   double get_Mhh(int i) const { return Mhh(i); }
   const Eigen::Array<double,2,1>& get_MAh() const { return MAh; }
   double get_MAh(int i) const { return MAh(i); }
   const Eigen::Array<double,2,1>& get_MHm() const { return MHm; }
   double get_MHm(int i) const { return MHm(i); }
   const Eigen::Array<double,3,1>& get_MFd() const { return MFd; }
   double get_MFd(int i) const { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const { return MFu; }
   double get_MFu(int i) const { return MFu(i); }
   const Eigen::Array<double,3,1>& get_MFe() const { return MFe; }
   double get_MFe(int i) const { return MFe(i); }
   double get_MVWm() const { return MVWm; }
   double get_MVP() const { return MVP; }
   double get_MVZ() const { return MVZ; }

   Eigen::Array<double,1,1> get_MChargedHiggs() const;
   Eigen::Array<double,1,1> get_MPseudoscalarHiggs() const;

   const Eigen::Matrix<double,2,2>& get_ZH() const { return ZH; }
   double get_ZH(int i, int k) const { return ZH(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZA() const { return ZA; }
   double get_ZA(int i, int k) const { return ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP() const { return ZP; }
   double get_ZP(int i, int k) const { return ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd() const { return Vd; }
   std::complex<double> get_Vd(int i, int k) const { return Vd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud() const { return Ud; }
   std::complex<double> get_Ud(int i, int k) const { return Ud(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu() const { return Vu; }
   std::complex<double> get_Vu(int i, int k) const { return Vu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu() const { return Uu; }
   std::complex<double> get_Uu(int i, int k) const { return Uu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve() const { return Ve; }
   std::complex<double> get_Ve(int i, int k) const { return Ve(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue() const { return Ue; }
   std::complex<double> get_Ue(int i, int k) const { return Ue(i,k); }

   double get_mass_matrix_VG() const;
   void calculate_MVG();
   double get_mass_matrix_VP() const;
   void calculate_MVP();
   double get_mass_matrix_VZ() const;
   void calculate_MVZ();
   double get_mass_matrix_VWm() const;
   void calculate_MVWm();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const;
   void calculate_MFv();
   Eigen::Matrix<double,2,2> get_mass_matrix_hh() const;
   void calculate_Mhh();
   Eigen::Matrix<double,2,2> get_mass_matrix_Ah() const;
   void calculate_MAh();
   Eigen::Matrix<double,2,2> get_mass_matrix_Hm() const;
   void calculate_MHm();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const;
   void calculate_MFd();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const;
   void calculate_MFu();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const;
   void calculate_MFe();

   double get_ewsb_eq_hh_1() const;
   double get_ewsb_eq_hh_2() const;

   double get_tan_beta() const;
   double get_beta() const;
   double get_alpha_h() const;  ///< CP-even Higgs mixing angle
   double get_alpha_em() const; ///< electromagnetic coupling
   double get_eta() const;      ///< deviation of CP-even Higgs mixing angle from SM limit
   double get_LambdaFive() const; ///< capital Lambda5, Eq (14) arxiv:1607.06292
   double ThetaW() const;
   double get_v() const;
   double get_v_sqr() const;

private:
   bool force_output{false};        ///< switch to force output of pole masses
   GeneralTHDM_physical physical{}; ///< contains the pole masses and mixings
   GeneralTHDM_problems problems{}; ///< problems

   // DR-bar masses
   double MVG{0.0};
   double MVWm{0.0};
   double MVP{0.0};
   double MVZ{0.0};
   Eigen::Array<double,3,1> MFv{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,2,1> Mhh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MAh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MHm{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,3,1> MFd{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFu{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFe{Eigen::Array<double,3,1>::Zero()};

   // DR-bar mixing matrices
   Eigen::Matrix<double,2,2> ZH{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZA{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZP{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vd{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ud{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Uu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ve{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ue{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
};

std::ostream& operator<<(std::ostream&, const GeneralTHDM_mass_eigenstates&);

} // namespace gm2calc

#endif
