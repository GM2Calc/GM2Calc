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
 * @file GeneralTHDM_mass_eigenstates.cpp
 * @brief implementation of the GeneralTHDM model class
 *
 * Contains the definition of the GeneralTHDM model class methods
 * which solve EWSB and calculate masses and mixings from MSbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.3 .
 */

#include "gm2calc/GeneralTHDM_mass_eigenstates.hpp"

#include "gm2_eigen_utils.hpp"
#include "gm2_linalg.hpp"
#include "gm2_numerics.hpp"
#include "gm2_raii.hpp"

#include <cmath>
#include <iostream>

namespace gm2calc {

namespace {

const double pi = 3.1415926535897932;
const double sqrt_inv = 0.70710678118654752; // 1/Sqrt[2]
const double gut_normalization = 0.77459666924148338; // Sqrt[3/5]

} // anonymous namespace

#define CLASSNAME GeneralTHDM_mass_eigenstates

void CLASSNAME::do_force_output(bool flag)
{
   force_output = flag;
}

bool CLASSNAME::do_force_output() const
{
   return force_output;
}

const GeneralTHDM_problems& CLASSNAME::get_problems() const
{
   return problems;
}

GeneralTHDM_problems& CLASSNAME::get_problems()
{
   return problems;
}

int CLASSNAME::solve_ewsb_tree_level()
{
   int error = 0;

   const double old_M112 = M112;
   const double old_M222 = M222;

   M112 = (0.25*(2*M122*v2 + 2*v2*M122 - 2*Lambda1*cube(v1) - Lambda7*
      cube(v2) - Lambda7*cube(v2) - 3*Lambda6*v2*sqr(v1) - 3*v2*Lambda6*sqr(v1)
      - 2*Lambda3*v1*sqr(v2) - 2*Lambda4*v1*sqr(v2) - Lambda5*v1*sqr(v2)
      - v1*Lambda5*sqr(v2)))/v1;
   M222 = (0.25*(2*M122*v1 + 2*v1*M122 - Lambda6*cube(v1) - Lambda6
      *cube(v1) - 2*Lambda2*cube(v2) - 2*Lambda3*v2*sqr(v1) - 2*Lambda4*v2*sqr(v1)
      - Lambda5*v2*sqr(v1) - v2*Lambda5*sqr(v1) - 3*Lambda7*v1*sqr(v2) - 3*
      v1*Lambda7*sqr(v2)))/v2;

   const bool is_finite = std::isfinite(M112) && std::isfinite(M222);

   if (!is_finite) {
      M112 = old_M112;
      M222 = old_M222;
      error = 1;
   }

   return error;
}

int CLASSNAME::solve_ewsb()
{
   return solve_ewsb_tree_level();
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "GeneralTHDM\n"
           "========================================\n";
   GeneralTHDM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level MSbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHm = " << MHm.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';

   ostr << "----------------------------------------\n"
           "tree-level MSbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';

   ostr << "----------------------------------------\n"
           "derived parameters:\n"
           "----------------------------------------\n";
   ostr << "v = " << get_v() << '\n';
   ostr << "theta_w = " << ThetaW() << '\n';
   ostr << "alpha_h = " << get_alpha_h() << '\n';
   ostr << "beta = " << get_beta() << '\n';
   ostr << "eta = " << get_eta() << '\n';
   ostr << "tan(beta) = " << get_tan_beta() << '\n';
}

/**
 * routine which finds the MSbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_MSbar_masses()
{
   const auto save_M112_raii = make_raii_save(M112);
   const auto save_M222_raii = make_raii_save(M222);

   solve_ewsb_tree_level();

   calculate_MVZ();
   calculate_MVWm();
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_MHm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MFv();
}

/**
 * reorders MSbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_MSbar_masses()
{
   move_goldstone_to(0, MVZ, MAh, ZA);
   move_goldstone_to(0, MVWm, MHm, ZP);

}

Eigen::Array<double,1,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHm_goldstone;
   MHm_goldstone(0) = MVWm;

   return remove_if_equal<double,2,1>(MHm, MHm_goldstone);
}

Eigen::Array<double,1,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,1,1> MAh_goldstone;
   MAh_goldstone(0) = MVZ;

   return remove_if_equal<double,2,1>(MAh, MAh_goldstone);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fv() const
{
   return Eigen::Matrix<double,3,3>::Zero();
}

void CLASSNAME::calculate_MFv()
{
   MFv.setZero();
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_hh() const
{
   Eigen::Matrix<double,2,2> mass_matrix_hh;

   mass_matrix_hh(0,0) = M112 + 1.5*Lambda6*v1*v2 + 1.5*v1*v2*Lambda6 +
      1.5*Lambda1*sqr(v1) + 0.5*Lambda3*sqr(v2) + 0.5*Lambda4*sqr(v2) + 0.25*
      Lambda5*sqr(v2) + 0.25*Lambda5*sqr(v2);
   mass_matrix_hh(0,1) = -0.5*M122 + Lambda3*v1*v2 + Lambda4*v1*v2 + 0.5*
      Lambda5*v1*v2 + 0.5*v1*v2*Lambda5 - 0.5*M122 + 0.75*Lambda6*
      sqr(v1) + 0.75*Lambda6*sqr(v1) + 0.75*Lambda7*sqr(v2) + 0.75*
      Lambda7*sqr(v2);
   mass_matrix_hh(1,1) = M222 + 1.5*Lambda7*v1*v2 + 1.5*v1*v2*Lambda7 +
      0.5*Lambda3*sqr(v1) + 0.5*Lambda4*sqr(v1) + 0.25*Lambda5*sqr(v1) + 0.25*
      Lambda5*sqr(v1) + 1.5*Lambda2*sqr(v2);

   symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());
   fs_diagonalize_hermitian<double,double,2>(mass_matrix_hh, Mhh, ZH);
   normalize_to_interval<2,2>(ZH);

   if (Mhh.minCoeff() < 0.) {
      problems.flag_tachyon("hh");
   }

   Mhh = sqrt(Mhh.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ah() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = M112 + 0.5*Lambda6*v1*v2 + 0.5*v1*v2*Lambda6 +
      0.5*Lambda1*sqr(v1) + 0.3872983346207417*g1*g2*std::cos(ThetaW())*std::sin(ThetaW()
      )*sqr(v1) + 0.5*Lambda3*sqr(v2) + 0.5*Lambda4*sqr(v2) - 0.25*Lambda5*sqr(
      v2) - 0.25*Lambda5*sqr(v2) + 0.25*sqr(g2)*sqr(v1)*sqr(std::cos(ThetaW())
      ) + 0.15*sqr(g1)*sqr(v1)*sqr(std::sin(ThetaW()));
   mass_matrix_Ah(0,1) = -0.5*M122 + 0.5*Lambda5*v1*v2 + 0.5*v1*v2*Lambda5
      - 0.5*M122 + 0.3872983346207417*g1*g2*v1*v2*std::cos(ThetaW())*std::sin(
      ThetaW()) + 0.25*Lambda6*sqr(v1) + 0.25*Lambda6*sqr(v1) + 0.25*
      Lambda7*sqr(v2) + 0.25*Lambda7*sqr(v2) + 0.25*v1*v2*sqr(g2)*sqr(std::cos
      (ThetaW())) + 0.15*v1*v2*sqr(g1)*sqr(std::sin(ThetaW()));
   mass_matrix_Ah(1,1) = M222 + 0.5*Lambda7*v1*v2 + 0.5*v1*v2*Lambda7 +
      0.5*Lambda3*sqr(v1) + 0.5*Lambda4*sqr(v1) - 0.25*Lambda5*sqr(v1) - 0.25*
      Lambda5*sqr(v1) + 0.5*Lambda2*sqr(v2) + 0.3872983346207417*g1*g2*
      std::cos(ThetaW())*std::sin(ThetaW())*sqr(v2) + 0.25*sqr(g2)*sqr(v2)*sqr(std::cos(ThetaW
      ())) + 0.15*sqr(g1)*sqr(v2)*sqr(std::sin(ThetaW()));

   symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());
   fs_diagonalize_hermitian<double,double,2>(mass_matrix_Ah, MAh, ZA);
   normalize_to_interval<2,2>(ZA);

   if (MAh.minCoeff() < 0.) {
      problems.flag_tachyon("Ah");
   }

   MAh = sqrt(MAh.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hm() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Hm;

   mass_matrix_Hm(0,0) = M112 + 0.5*Lambda6*v1*v2 + 0.5*v1*v2*Lambda6 +
      0.5*Lambda1*sqr(v1) + 0.25*sqr(g2)*sqr(v1) + 0.5*Lambda3*sqr(v2);
   mass_matrix_Hm(0,1) = 0.5*Lambda4*v1*v2 + 0.5*Lambda5*v1*v2 - M122 +
      0.25*v1*v2*sqr(g2) + 0.5*Lambda6*sqr(v1) + 0.5*Lambda7*sqr(v2
      );
   mass_matrix_Hm(1,0) = -M122 + 0.5*Lambda4*v1*v2 + 0.5*v1*v2*Lambda5 +
      0.25*v1*v2*sqr(g2) + 0.5*Lambda6*sqr(v1) + 0.5*Lambda7*sqr(v2);
   mass_matrix_Hm(1,1) = M222 + 0.5*Lambda7*v1*v2 + 0.5*v1*v2*Lambda7 +
      0.5*Lambda3*sqr(v1) + 0.5*Lambda2*sqr(v2) + 0.25*sqr(g2)*sqr(v2);

   return mass_matrix_Hm;
}

void CLASSNAME::calculate_MHm()
{
   const auto mass_matrix_Hm(get_mass_matrix_Hm());
   fs_diagonalize_hermitian<double,double,2>(mass_matrix_Hm, MHm, ZP);
   normalize_to_interval<2,2>(ZP);

   if (MHm.minCoeff() < 0.) {
      problems.flag_tachyon("Hm");
   }

   MHm = sqrt(MHm.cwiseAbs());
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{
   return sqrt_inv*(v1*Yd + v2*Xd);
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());
   fs_svd<double,3,3>(mass_matrix_Fd, MFd, Vd, Ud);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{
   return sqrt_inv*(v2*Yu + v1*Xu);
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());
   fs_svd<double,3,3>(mass_matrix_Fu, MFu, Vu, Uu);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{
   return sqrt_inv*(v1*Ye + v2*Xe);
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());
   fs_svd<double,3,3>(mass_matrix_Fe, MFe, Ve, Ue);
}

double CLASSNAME::get_mass_matrix_VWm() const
{
   const double mass_matrix_VWm = 0.25*sqr(g2)*(sqr(v1) + sqr(v2));
   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{
   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = abs_sqrt(mass_matrix_VWm);
}

double CLASSNAME::get_mass_matrix_VZ() const
{
   const double mass_matrix_VZ = 0.25*(sqr(v1) + sqr(v2))*sqr(g2*std::cos(
      ThetaW()) + gut_normalization*g1*std::sin(ThetaW()));

   return mass_matrix_VZ;
}

void CLASSNAME::calculate_MVZ()
{
   const auto mass_matrix_VZ = get_mass_matrix_VZ();
   MVZ = abs_sqrt(mass_matrix_VZ);
}

double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = M112*v1 - 0.5*M122*v2 - 0.5*v2*M122 + 0.5*Lambda1*cube
      (v1) + 0.25*Lambda7*cube(v2) + 0.25*Lambda7*cube(v2) + 0.75*Lambda6*v2
      *sqr(v1) + 0.75*v2*Lambda6*sqr(v1) + 0.5*Lambda3*v1*sqr(v2) + 0.5*
      Lambda4*v1*sqr(v2) + 0.25*Lambda5*v1*sqr(v2) + 0.25*v1*Lambda5*sqr(v2);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   double result = -0.5*M122*v1 + M222*v2 - 0.5*v1*M122 + 0.25*Lambda6*
      cube(v1) + 0.25*Lambda6*cube(v1) + 0.5*Lambda2*cube(v2) + 0.5*Lambda3*
      v2*sqr(v1) + 0.5*Lambda4*v2*sqr(v1) + 0.25*Lambda5*v2*sqr(v1) + 0.25*v2*
      Lambda5*sqr(v1) + 0.75*Lambda7*v1*sqr(v2) + 0.75*v1*Lambda7*sqr(v2);

   return result;
}

/**
 * Returns the CP-odd Higgs mixing angle \f$\beta\f$.
 *
 * The mixing angle \f$\beta\f$ is chosen such that
 * \f$0 \leq \beta \leq \pi/2\f$.
 *
 * @return CP-odd Higgs mixing angle \f$\beta\f$
 */
double CLASSNAME::get_beta() const
{
   return std::atan(get_tan_beta());
}

double CLASSNAME::get_tan_beta() const
{
   return v2/v1;
}

/**
 * Returns the CP-even Higgs mixing angle \f$\alpha_h\f$.
 *
 * The mixing angle \f$\alpha_h\f$ is chosen such that
 * \f$-\pi/2 \leq \beta - \alpha_h \leq \pi/2\f$.
 *
 * @return CP-even Higgs mixing angle \f$\alpha_h\f$
 */
double CLASSNAME::get_alpha_h() const
{
   double alpha_h = std::atan2(ZH(1,1), ZH(0,1));
   const double bma = get_beta() - alpha_h;

   if (bma < -pi/2) { alpha_h -= pi; }
   if (bma >  pi/2) { alpha_h += pi; }

   return alpha_h;
}

double CLASSNAME::get_alpha_em() const
{
   const double gY = g1*gut_normalization;
   const double e2 = sqr(gY)*sqr(g2)/(sqr(gY) + sqr(g2));
   return e2/12.566370614359173; // e^2/(4 Pi)
}

/**
 * Returns the CP-even Higgs mixing angle contribution \f$\eta\f$ that
 * describes the deviation from the SM limit.
 *
 * The \f$\eta\f$ parameter is defined as
 * \f$\beta - \alpha = \pi/2 - \eta\f$.
 *
 * @return \f$\eta\f$
 */
double CLASSNAME::get_eta() const
{
   return pi/2 + get_alpha_h() - get_beta();
}

/**
 * Returns \f$\Lambda_5\f$, Eq (14) arxiv:1607.06292
 * @return \f$\Lambda_5\f$
 */
double CLASSNAME::get_LambdaFive() const
{
   const double tb = get_tan_beta();
   const double sbcb = tb/(1.0 + sqr(tb)); // Sin[beta] Cos[beta]
   return 2*get_M122()/(get_v_sqr() * sbcb);
}

double CLASSNAME::ThetaW() const
{
   return std::atan(gut_normalization*g1/g2);
}

double CLASSNAME::get_v() const
{
   return std::sqrt(get_v_sqr());
}

double CLASSNAME::get_v_sqr() const
{
   return sqr(v1) + sqr(v2);
}

void CLASSNAME::set_alpha_em_and_cw(double alpha_em, double cw)
{
   const double e = std::sqrt(alpha_em*4*pi);
   const double sw = std::sqrt(1 - sqr(cw));
   g1 = e/(cw*gut_normalization);
   g2 = e/sw;
}

void CLASSNAME::set_tan_beta_and_v(double tan_beta, double v)
{
   const double sq = std::sqrt(1.0 + sqr(tan_beta));
   const double sb = tan_beta/sq;
   const double cb = 1/sq;
   v1 = v*cb;
   v2 = v*sb;
}

std::ostream& operator<<(std::ostream& ostr, const CLASSNAME& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace gm2calc
