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
 * @file MSSMNoFV_onshell_mass_eigenstates.cpp
 * @brief implementation of the MSSMNoFV_onshell model class
 *
 * Contains the definition of the MSSMNoFV_onshell model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Wed 22 Jul 2015 18:14:31 with FlexibleSUSY
 * 1.2.1 (git commit: v1.2.1-49-gfc4c300) and SARAH 4.5.8 .
 */

#include "MSSMNoFV_onshell_mass_eigenstates.hpp"
#include "eigen_utils.hpp"
#include "linalg2.hpp"
#include "numerics2.hpp"
#include "ffunctions.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

namespace {
template <typename Derived>
void Symmetrize(Eigen::MatrixBase<Derived>& m)
{
   static_assert(Eigen::MatrixBase<Derived>::RowsAtCompileTime ==
                 Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                 "Symmetrize is only defined for squared matrices");

   for (int i = 0; i < Eigen::MatrixBase<Derived>::RowsAtCompileTime; i++)
      for (int k = 0; k < i; k++)
         m(i,k) = m(k,i);
}

template <typename Derived>
void Hermitianize(Eigen::MatrixBase<Derived>& m)
{
   static_assert(Eigen::MatrixBase<Derived>::RowsAtCompileTime ==
                 Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                 "Hermitianize is only defined for squared matrices");

   for (int i = 0; i < Eigen::MatrixBase<Derived>::RowsAtCompileTime; i++)
      for (int k = 0; k < i; k++)
         m(i,k) = std::conj(m(k,i));
}
} // anonymous namespace

namespace gm2calc {

using namespace flexiblesusy;

#define CLASSNAME MSSMNoFV_onshell_mass_eigenstates
#define PHYSICAL(parameter) physical.parameter
#define MODELPARAMETER(parameter) model->get_##parameter()

CLASSNAME::MSSMNoFV_onshell_mass_eigenstates()
   : MSSMNoFV_onshell_soft_parameters()
   , force_output(false)
   , precision(1.0e-3)
   , physical()
   , problems()
   , MVG(0), MGlu(0), MVP(0), MVZ(0), MFd(0), MFs(0), MFb(0), MFu(0), MFc(0),
      MFt(0), MFve(0), MFvm(0), MFvt(0), MFe(0), MFm(0), MFtau(0), MSveL(0), MSvmL
      (0), MSvtL(0), MSd(Eigen::Array<double,2,1>::Zero()), MSu(Eigen::Array<
      double,2,1>::Zero()), MSe(Eigen::Array<double,2,1>::Zero()), MSm(
      Eigen::Array<double,2,1>::Zero()), MStau(Eigen::Array<double,2,1>::Zero()),
      MSs(Eigen::Array<double,2,1>::Zero()), MSc(Eigen::Array<double,2,1>::Zero())
      , MSb(Eigen::Array<double,2,1>::Zero()), MSt(Eigen::Array<double,2,1>::Zero(
      )), Mhh(Eigen::Array<double,2,1>::Zero()), MAh(Eigen::Array<double,2,1>
      ::Zero()), MHpm(Eigen::Array<double,2,1>::Zero()), MChi(Eigen::Array<double,
      4,1>::Zero()), MCha(Eigen::Array<double,2,1>::Zero()), MVWm(0)

   , ZD(Eigen::Matrix<double,2,2>::Zero()), ZU(Eigen::Matrix<double,2,2>::Zero(
      )), ZE(Eigen::Matrix<double,2,2>::Zero()), ZM(Eigen::Matrix<double,2,2>
      ::Zero()), ZTau(Eigen::Matrix<double,2,2>::Zero()), ZS(Eigen::Matrix<double,
      2,2>::Zero()), ZC(Eigen::Matrix<double,2,2>::Zero()), ZB(Eigen::Matrix<
      double,2,2>::Zero()), ZT(Eigen::Matrix<double,2,2>::Zero()), ZH(
      Eigen::Matrix<double,2,2>::Zero()), ZA(Eigen::Matrix<double,2,2>::Zero()),
      ZP(Eigen::Matrix<double,2,2>::Zero()), ZN(Eigen::Matrix<std::complex<double>
      ,4,4>::Zero()), UM(Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP(
      Eigen::Matrix<std::complex<double>,2,2>::Zero())

   , PhaseGlu(1,0)

{
}

CLASSNAME::~MSSMNoFV_onshell_mass_eigenstates()
{
}

void CLASSNAME::do_force_output(bool flag)
{
   force_output = flag;
}

bool CLASSNAME::do_force_output() const
{
   return force_output;
}

const MSSMNoFV_onshell_physical& CLASSNAME::get_physical() const
{
   return physical;
}

MSSMNoFV_onshell_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const MSSMNoFV_onshell_physical& physical_)
{
   physical = physical_;
}

const MSSMNoFV_onshell_problems& CLASSNAME::get_problems() const
{
   return problems;
}

MSSMNoFV_onshell_problems& CLASSNAME::get_problems()
{
   return problems;
}

int CLASSNAME::solve_ewsb_tree_level()
{
   return solve_ewsb_tree_level_via_soft_higgs_masses();
}

int CLASSNAME::solve_ewsb_tree_level_via_soft_higgs_masses()
{
   int error = 0;

   const double new_mHd2 = (0.025*(-40*vd*sqr(Mu) + 20*vu*BMu + 20*vu*
      BMu - 3*std::pow(vd,3)*sqr(g1) - 5*std::pow(vd,3)*sqr(g2) + 3*vd*sqr(g1)*sqr
      (vu) + 5*vd*sqr(g2)*sqr(vu)))/vd;
   const double new_mHu2 = (0.025*(-40*vu*sqr(Mu) + 20*vd*BMu + 20*vd*
      BMu - 3*std::pow(vu,3)*sqr(g1) - 5*std::pow(vu,3)*sqr(g2) + 3*vu*sqr(g1)*sqr
      (vd) + 5*vu*sqr(g2)*sqr(vd)))/vu;

   if (std::isfinite(new_mHd2))
      mHd2 = new_mHd2;
   else
      error = 1;

   if (std::isfinite(new_mHu2))
      mHu2 = new_mHu2;
   else
      error = 1;


   return error;
}

int CLASSNAME::solve_ewsb()
{
   return solve_ewsb_tree_level();
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "MSSMNoFV_onshell\n"
           "========================================\n";
   MSSMNoFV_onshell_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MFd = " << MFd << '\n';
   ostr << "MFs = " << MFs << '\n';
   ostr << "MFb = " << MFb << '\n';
   ostr << "MFu = " << MFu << '\n';
   ostr << "MFc = " << MFc << '\n';
   ostr << "MFt = " << MFt << '\n';
   ostr << "MFve = " << MFve << '\n';
   ostr << "MFvm = " << MFvm << '\n';
   ostr << "MFvt = " << MFvt << '\n';
   ostr << "MFe = " << MFe << '\n';
   ostr << "MFm = " << MFm << '\n';
   ostr << "MFtau = " << MFtau << '\n';
   ostr << "MSveL = " << MSveL << '\n';
   ostr << "MSvmL = " << MSvmL << '\n';
   ostr << "MSvtL = " << MSvtL << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "MSm = " << MSm.transpose() << '\n';
   ostr << "MStau = " << MStau.transpose() << '\n';
   ostr << "MSs = " << MSs.transpose() << '\n';
   ostr << "MSc = " << MSc.transpose() << '\n';
   ostr << "MSb = " << MSb.transpose() << '\n';
   ostr << "MSt = " << MSt.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';

   ostr << "----------------------------------------\n"
           "tree-level DRbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZD = " << ZD << '\n';
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZM = " << ZM << '\n';
   ostr << "ZTau = " << ZTau << '\n';
   ostr << "ZS = " << ZS << '\n';
   ostr << "ZC = " << ZC << '\n';
   ostr << "ZB = " << ZB << '\n';
   ostr << "ZT = " << ZT << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';

   physical.print(ostr);
   problems.print(ostr);
}

/**
 * routine which finds the DRbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{
   const auto old_mHd2 = mHd2;
   const auto old_mHu2 = mHu2;

   solve_ewsb_tree_level_via_soft_higgs_masses();

   calculate_MVG();
   calculate_MVP();
   calculate_MVZ();
   calculate_MVWm();
   calculate_MGlu();
   calculate_MFd();
   calculate_MFs();
   calculate_MFb();
   calculate_MFu();
   calculate_MFc();
   calculate_MFt();
   calculate_MFve();
   calculate_MFvm();
   calculate_MFvt();
   calculate_MFe();
   calculate_MFm();
   calculate_MFtau();
   calculate_MSveL();
   calculate_MSvmL();
   calculate_MSvtL();
   calculate_MSd();
   calculate_MSu();
   calculate_MSe();
   calculate_MSm();
   calculate_MStau();
   calculate_MSs();
   calculate_MSc();
   calculate_MSb();
   calculate_MSt();
   calculate_Mhh();
   calculate_MAh();
   calculate_MHpm();
   calculate_MChi();
   calculate_MCha();

   mHd2 = old_mHd2;
   mHu2 = old_mHu2;

}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MGlu) = MGlu;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(MFs) = MFs;
   PHYSICAL(MFb) = MFb;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(MFc) = MFc;
   PHYSICAL(MFt) = MFt;
   PHYSICAL(MFve) = MFve;
   PHYSICAL(MFvm) = MFvm;
   PHYSICAL(MFvt) = MFvt;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(MFm) = MFm;
   PHYSICAL(MFtau) = MFtau;
   PHYSICAL(MSveL) = MSveL;
   PHYSICAL(MSvmL) = MSvmL;
   PHYSICAL(MSvtL) = MSvtL;
   PHYSICAL(MSd) = MSd;
   PHYSICAL(ZD) = ZD;
   PHYSICAL(MSu) = MSu;
   PHYSICAL(ZU) = ZU;
   PHYSICAL(MSe) = MSe;
   PHYSICAL(ZE) = ZE;
   PHYSICAL(MSm) = MSm;
   PHYSICAL(ZM) = ZM;
   PHYSICAL(MStau) = MStau;
   PHYSICAL(ZTau) = ZTau;
   PHYSICAL(MSs) = MSs;
   PHYSICAL(ZS) = ZS;
   PHYSICAL(MSc) = MSc;
   PHYSICAL(ZC) = ZC;
   PHYSICAL(MSb) = MSb;
   PHYSICAL(ZB) = ZB;
   PHYSICAL(MSt) = MSt;
   PHYSICAL(ZT) = ZT;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(ZH) = ZH;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(ZA) = ZA;
   PHYSICAL(MHpm) = MHpm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MChi) = MChi;
   PHYSICAL(ZN) = ZN;
   PHYSICAL(MCha) = MCha;
   PHYSICAL(UM) = UM;
   PHYSICAL(UP) = UP;
   PHYSICAL(MVWm) = MVWm;

}

/**
 * reorders DRbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_DRbar_masses()
{
   move_goldstone_to(0, MVZ, MAh, ZA);
   move_goldstone_to(0, MVWm, MHpm, ZP);
}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{
   move_goldstone_to(0, MVZ, PHYSICAL(MAh), PHYSICAL(ZA));
   move_goldstone_to(0, MVWm, PHYSICAL(MHpm), PHYSICAL(ZP));
}

void CLASSNAME::clear_DRbar_parameters()
{
   MVG = 0.;
   MGlu = 0.;
   MVP = 0.;
   MVZ = 0.;
   MFd = 0.;
   MFs = 0.;
   MFb = 0.;
   MFu = 0.;
   MFc = 0.;
   MFt = 0.;
   MFve = 0.;
   MFvm = 0.;
   MFvt = 0.;
   MFe = 0.;
   MFm = 0.;
   MFtau = 0.;
   MSveL = 0.;
   MSvmL = 0.;
   MSvtL = 0.;
   MSd = Eigen::Matrix<double,2,1>::Zero();
   ZD = Eigen::Matrix<double,2,2>::Zero();
   MSu = Eigen::Matrix<double,2,1>::Zero();
   ZU = Eigen::Matrix<double,2,2>::Zero();
   MSe = Eigen::Matrix<double,2,1>::Zero();
   ZE = Eigen::Matrix<double,2,2>::Zero();
   MSm = Eigen::Matrix<double,2,1>::Zero();
   ZM = Eigen::Matrix<double,2,2>::Zero();
   MStau = Eigen::Matrix<double,2,1>::Zero();
   ZTau = Eigen::Matrix<double,2,2>::Zero();
   MSs = Eigen::Matrix<double,2,1>::Zero();
   ZS = Eigen::Matrix<double,2,2>::Zero();
   MSc = Eigen::Matrix<double,2,1>::Zero();
   ZC = Eigen::Matrix<double,2,2>::Zero();
   MSb = Eigen::Matrix<double,2,1>::Zero();
   ZB = Eigen::Matrix<double,2,2>::Zero();
   MSt = Eigen::Matrix<double,2,1>::Zero();
   ZT = Eigen::Matrix<double,2,2>::Zero();
   Mhh = Eigen::Matrix<double,2,1>::Zero();
   ZH = Eigen::Matrix<double,2,2>::Zero();
   MAh = Eigen::Matrix<double,2,1>::Zero();
   ZA = Eigen::Matrix<double,2,2>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,4,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MCha = Eigen::Matrix<double,2,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MVWm = 0.;

   PhaseGlu = std::complex<double>(1.,0.);

}

void CLASSNAME::clear()
{
   MSSMNoFV_onshell_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

std::string CLASSNAME::name() const
{
   return "MSSMNoFV_onshell";
}

Eigen::Array<double,1,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHpm_ChargedHiggs;
   Eigen::Array<double,1,1> MHpm_goldstone;

   MHpm_goldstone(0) = MVWm;

   remove_if_equal(MHpm, MHpm_goldstone, MHpm_ChargedHiggs);

   return MHpm_ChargedHiggs;
}

Eigen::Array<double,1,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,1,1> MAh_PseudoscalarHiggs;
   Eigen::Array<double,1,1> MAh_goldstone;

   MAh_goldstone(0) = MVZ;

   remove_if_equal(MAh, MAh_goldstone, MAh_PseudoscalarHiggs);

   return MAh_PseudoscalarHiggs;
}




double CLASSNAME::get_mass_matrix_VG() const
{
   return 0;
}

void CLASSNAME::calculate_MVG()
{
   MVG = 0;
}

double CLASSNAME::get_mass_matrix_Glu() const
{
   return MassG;
}

void CLASSNAME::calculate_MGlu()
{
   const auto mass_matrix_Glu = get_mass_matrix_Glu();
   PhaseGlu = std::polar(1., 0.5 * std::arg(std::complex<double>(mass_matrix_Glu)));
   MGlu = std::abs(mass_matrix_Glu);
}

double CLASSNAME::get_mass_matrix_VP() const
{
   return 0;
}

void CLASSNAME::calculate_MVP()
{
   const auto mass_matrix_VP = get_mass_matrix_VP();
   MVP = std::abs(mass_matrix_VP);
}

double CLASSNAME::get_mass_matrix_VZ() const
{
   const double mass_matrix_VZ = 0.25*(sqr(vd) + sqr(vu))*sqr(g2*std::cos(
      ThetaW()) + 0.7745966692414834*g1*std::sin(ThetaW()));

   return mass_matrix_VZ;
}

void CLASSNAME::calculate_MVZ()
{
   const auto mass_matrix_VZ = get_mass_matrix_VZ();
   MVZ = std::abs(mass_matrix_VZ);

   MVZ = std::sqrt(std::abs(MVZ));
}

double CLASSNAME::get_mass_matrix_Fd() const
{
   const double mass_matrix_Fd = 0.7071067811865475*vd*Yd(0,0);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd = get_mass_matrix_Fd();
   MFd = std::abs(mass_matrix_Fd);
}

double CLASSNAME::get_mass_matrix_Fs() const
{
   const double mass_matrix_Fs = 0.7071067811865475*vd*Yd(1,1);

   return mass_matrix_Fs;
}

void CLASSNAME::calculate_MFs()
{
   const auto mass_matrix_Fs = get_mass_matrix_Fs();
   MFs = std::abs(mass_matrix_Fs);
}

double CLASSNAME::get_mass_matrix_Fb() const
{
   const double mass_matrix_Fb = 0.7071067811865475*vd*Yd(2,2);

   return mass_matrix_Fb;
}

void CLASSNAME::calculate_MFb()
{
   const auto mass_matrix_Fb = get_mass_matrix_Fb();
   MFb = std::abs(mass_matrix_Fb);
}

double CLASSNAME::get_mass_matrix_Fu() const
{
   const double mass_matrix_Fu = 0.7071067811865475*vu*Yu(0,0);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu = get_mass_matrix_Fu();
   MFu = std::abs(mass_matrix_Fu);
}

double CLASSNAME::get_mass_matrix_Fc() const
{
   const double mass_matrix_Fc = 0.7071067811865475*vu*Yu(1,1);

   return mass_matrix_Fc;
}

void CLASSNAME::calculate_MFc()
{
   const auto mass_matrix_Fc = get_mass_matrix_Fc();
   MFc = std::abs(mass_matrix_Fc);
}

double CLASSNAME::get_mass_matrix_Ft() const
{
   const double mass_matrix_Ft = 0.7071067811865475*vu*Yu(2,2);

   return mass_matrix_Ft;
}

void CLASSNAME::calculate_MFt()
{
   const auto mass_matrix_Ft = get_mass_matrix_Ft();
   MFt = std::abs(mass_matrix_Ft);
}

double CLASSNAME::get_mass_matrix_Fve() const
{
   return 0;
}

void CLASSNAME::calculate_MFve()
{
   const auto mass_matrix_Fve = get_mass_matrix_Fve();
   MFve = std::abs(mass_matrix_Fve);
}

double CLASSNAME::get_mass_matrix_Fvm() const
{
   return 0;
}

void CLASSNAME::calculate_MFvm()
{
   const auto mass_matrix_Fvm = get_mass_matrix_Fvm();
   MFvm = std::abs(mass_matrix_Fvm);
}

double CLASSNAME::get_mass_matrix_Fvt() const
{
   return 0;
}

void CLASSNAME::calculate_MFvt()
{
   const auto mass_matrix_Fvt = get_mass_matrix_Fvt();
   MFvt = std::abs(mass_matrix_Fvt);
}

double CLASSNAME::get_mass_matrix_Fe() const
{
   const double mass_matrix_Fe = 0.7071067811865475*vd*Ye(0,0);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe = get_mass_matrix_Fe();
   MFe = std::abs(mass_matrix_Fe);
}

double CLASSNAME::get_mass_matrix_Fm() const
{
   const double mass_matrix_Fm = 0.7071067811865475*vd*Ye(1,1);

   return mass_matrix_Fm;
}

void CLASSNAME::calculate_MFm()
{
   const auto mass_matrix_Fm = get_mass_matrix_Fm();
   MFm = std::abs(mass_matrix_Fm);
}

double CLASSNAME::get_mass_matrix_Ftau() const
{
   const double mass_matrix_Ftau = 0.7071067811865475*vd*Ye(2,2);

   return mass_matrix_Ftau;
}

void CLASSNAME::calculate_MFtau()
{
   const auto mass_matrix_Ftau = get_mass_matrix_Ftau();
   MFtau = std::abs(mass_matrix_Ftau);
}

double CLASSNAME::get_mass_matrix_SveL() const
{
   const double mass_matrix_SveL = 0.125*(8*ml2(0,0) - 0.6*sqr(g1)*(
      -sqr(vd) + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));

   return mass_matrix_SveL;
}

void CLASSNAME::calculate_MSveL()
{
   const auto mass_matrix_SveL = get_mass_matrix_SveL();
   MSveL = std::abs(mass_matrix_SveL);

   // if (MSveL < 0.)
   //    problems.flag_tachyon("SveL");

   MSveL = std::sqrt(std::abs(MSveL));
}

double CLASSNAME::get_mass_matrix_SvmL() const
{
   const double mass_matrix_SvmL = 0.125*(8*ml2(1,1) - 0.6*sqr(g1)*(
      -sqr(vd) + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));

   return mass_matrix_SvmL;
}

void CLASSNAME::calculate_MSvmL()
{
   const auto mass_matrix_SvmL = get_mass_matrix_SvmL();
   MSvmL = std::abs(mass_matrix_SvmL);

   if (MSvmL < 0.)
      problems.flag_tachyon("SvmL");

   MSvmL = std::sqrt(std::abs(MSvmL));
}

double CLASSNAME::get_mass_matrix_SvtL() const
{
   const double mass_matrix_SvtL = 0.125*(8*ml2(2,2) - 0.6*sqr(g1)*(
      -sqr(vd) + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));

   return mass_matrix_SvtL;
}

void CLASSNAME::calculate_MSvtL()
{
   const auto mass_matrix_SvtL = get_mass_matrix_SvtL();
   MSvtL = std::abs(mass_matrix_SvtL);

   // if (MSvtL < 0.)
   //    problems.flag_tachyon("SvtL");

   MSvtL = std::sqrt(std::abs(MSvtL));
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sd() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*sqr(Yd(0,0))*sqr(vd) - 0.025*
      sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) + 0.025*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Sd(0,1) = 0.7071067811865475*vd*std::conj(TYd(0,0)) -
      0.7071067811865475*vu*std::conj(Yd(0,0))*Mu;
   mass_matrix_Sd(1,1) = md2(0,0) + 0.5*sqr(Yd(0,0))*sqr(vd) - 0.05*
      sqr(g1)*sqr(vd) + 0.05*sqr(g1)*sqr(vu);

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);

   // if (MSd.minCoeff() < 0.)
   //    problems.flag_tachyon("Sd");

   MSd = sqrt(MSd.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Su() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*sqr(g1)*sqr(vd) + 0.125*sqr(g2)
      *sqr(vd) + 0.5*sqr(Yu(0,0))*sqr(vu) + 0.025*sqr(g1)*sqr(vu) - 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Su(0,1) = 0.7071067811865475*vu*std::conj(TYu(0,0)) -
      0.7071067811865475*vd*std::conj(Yu(0,0))*Mu;
   mass_matrix_Su(1,1) = mu2(0,0) + 0.1*sqr(g1)*sqr(vd) + 0.5*sqr(Yu(0
      ,0))*sqr(vu) - 0.1*sqr(g1)*sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);

   // if (MSu.minCoeff() < 0.)
   //    problems.flag_tachyon("Su");

   MSu = sqrt(MSu.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Se() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*sqr(Ye(0,0))*sqr(vd) + 0.075*
      sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Se(0,1) = 0.7071067811865475*vd*std::conj(TYe(0,0)) -
      0.7071067811865475*vu*std::conj(Ye(0,0))*Mu;
   mass_matrix_Se(1,1) = me2(0,0) + 0.5*sqr(Ye(0,0))*sqr(vd) - 0.15*
      sqr(g1)*sqr(vd) + 0.15*sqr(g1)*sqr(vu);

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);

   // if (MSe.minCoeff() < 0.)
   //    problems.flag_tachyon("Se");

   MSe = sqrt(MSe.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sm() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Sm;

   mass_matrix_Sm(0,0) = ml2(1,1) + 0.5*sqr(Ye(1,1))*sqr(vd) + 0.075*
      sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Sm(0,1) = 0.7071067811865475*vd*std::conj(TYe(1,1)) -
      0.7071067811865475*vu*std::conj(Ye(1,1))*Mu;
   mass_matrix_Sm(1,1) = me2(1,1) + 0.5*sqr(Ye(1,1))*sqr(vd) - 0.15*
      sqr(g1)*sqr(vd) + 0.15*sqr(g1)*sqr(vu);

   Hermitianize(mass_matrix_Sm);

   return mass_matrix_Sm;
}

void CLASSNAME::calculate_MSm()
{
   const auto mass_matrix_Sm(get_mass_matrix_Sm());
   fs_diagonalize_hermitian(mass_matrix_Sm, MSm, ZM);

   if (MSm.minCoeff() < 0.)
      problems.flag_tachyon("Sm");

   MSm = sqrt(MSm.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Stau() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Stau;

   mass_matrix_Stau(0,0) = ml2(2,2) + 0.5*sqr(Ye(2,2))*sqr(vd) + 0.075
      *sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Stau(0,1) = 0.7071067811865475*vd*std::conj(TYe(2,2)) -
      0.7071067811865475*vu*std::conj(Ye(2,2))*Mu;
   mass_matrix_Stau(1,1) = me2(2,2) + 0.5*sqr(Ye(2,2))*sqr(vd) - 0.15*
      sqr(g1)*sqr(vd) + 0.15*sqr(g1)*sqr(vu);

   Hermitianize(mass_matrix_Stau);

   return mass_matrix_Stau;
}

void CLASSNAME::calculate_MStau()
{
   const auto mass_matrix_Stau(get_mass_matrix_Stau());
   fs_diagonalize_hermitian(mass_matrix_Stau, MStau, ZTau);

   if (MStau.minCoeff() < 0.)
      problems.flag_tachyon("Stau");

   MStau = sqrt(MStau.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ss() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Ss;

   mass_matrix_Ss(0,0) = mq2(1,1) + 0.5*sqr(Yd(1,1))*sqr(vd) - 0.025*
      sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) + 0.025*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Ss(0,1) = 0.7071067811865475*vd*std::conj(TYd(1,1)) -
      0.7071067811865475*vu*std::conj(Yd(1,1))*Mu;
   mass_matrix_Ss(1,1) = md2(1,1) + 0.5*sqr(Yd(1,1))*sqr(vd) - 0.05*
      sqr(g1)*sqr(vd) + 0.05*sqr(g1)*sqr(vu);

   Hermitianize(mass_matrix_Ss);

   return mass_matrix_Ss;
}

void CLASSNAME::calculate_MSs()
{
   const auto mass_matrix_Ss(get_mass_matrix_Ss());
   fs_diagonalize_hermitian(mass_matrix_Ss, MSs, ZS);

   // if (MSs.minCoeff() < 0.)
   //    problems.flag_tachyon("Ss");

   MSs = sqrt(MSs.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sc() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Sc;

   mass_matrix_Sc(0,0) = mq2(1,1) - 0.025*sqr(g1)*sqr(vd) + 0.125*sqr(g2)
      *sqr(vd) + 0.5*sqr(Yu(1,1))*sqr(vu) + 0.025*sqr(g1)*sqr(vu) - 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Sc(0,1) = 0.7071067811865475*vu*std::conj(TYu(1,1)) -
      0.7071067811865475*vd*std::conj(Yu(1,1))*Mu;
   mass_matrix_Sc(1,1) = mu2(1,1) + 0.1*sqr(g1)*sqr(vd) + 0.5*sqr(Yu(1
      ,1))*sqr(vu) - 0.1*sqr(g1)*sqr(vu);

   Hermitianize(mass_matrix_Sc);

   return mass_matrix_Sc;
}

void CLASSNAME::calculate_MSc()
{
   const auto mass_matrix_Sc(get_mass_matrix_Sc());
   fs_diagonalize_hermitian(mass_matrix_Sc, MSc, ZC);

   // if (MSc.minCoeff() < 0.)
   //    problems.flag_tachyon("Sc");

   MSc = sqrt(MSc.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sb() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Sb;

   mass_matrix_Sb(0,0) = mq2(2,2) + 0.5*sqr(Yd(2,2))*sqr(vd) - 0.025*
      sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) + 0.025*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Sb(0,1) = 0.7071067811865475*vd*std::conj(TYd(2,2)) -
      0.7071067811865475*vu*std::conj(Yd(2,2))*Mu;
   mass_matrix_Sb(1,1) = md2(2,2) + 0.5*sqr(Yd(2,2))*sqr(vd) - 0.05*
      sqr(g1)*sqr(vd) + 0.05*sqr(g1)*sqr(vu);

   Hermitianize(mass_matrix_Sb);

   return mass_matrix_Sb;
}

void CLASSNAME::calculate_MSb()
{
   const auto mass_matrix_Sb(get_mass_matrix_Sb());
   fs_diagonalize_hermitian(mass_matrix_Sb, MSb, ZB);

   if (MSb.minCoeff() < 0.)
      problems.flag_tachyon("Sb");

   MSb = sqrt(MSb.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_St() const
{
   Eigen::Matrix<double,2,2> mass_matrix_St;

   mass_matrix_St(0,0) = mq2(2,2) - 0.025*sqr(g1)*sqr(vd) + 0.125*sqr(g2)
      *sqr(vd) + 0.5*sqr(Yu(2,2))*sqr(vu) + 0.025*sqr(g1)*sqr(vu) - 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_St(0,1) = 0.7071067811865475*vu*std::conj(TYu(2,2)) -
      0.7071067811865475*vd*std::conj(Yu(2,2))*Mu;
   mass_matrix_St(1,1) = mu2(2,2) + 0.1*sqr(g1)*sqr(vd) + 0.5*sqr(Yu(2
      ,2))*sqr(vu) - 0.1*sqr(g1)*sqr(vu);

   Hermitianize(mass_matrix_St);

   return mass_matrix_St;
}

void CLASSNAME::calculate_MSt()
{
   const auto mass_matrix_St(get_mass_matrix_St());
   fs_diagonalize_hermitian(mass_matrix_St, MSt, ZT);

   if (MSt.minCoeff() < 0.)
      problems.flag_tachyon("St");

   MSt = sqrt(MSt.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_hh() const
{
   Eigen::Matrix<double,2,2> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + sqr(Mu) + 0.225*sqr(g1)*sqr(vd) +
      0.375*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) - 0.125*sqr(g2)*sqr(vu);
   mass_matrix_hh(0,1) = -0.5*BMu - 0.5*BMu - 0.15*vd*vu*sqr(g1) -
      0.25*vd*vu*sqr(g2);
   mass_matrix_hh(1,1) = mHu2 + sqr(Mu) - 0.075*sqr(g1)*sqr(vd) -
      0.125*sqr(g2)*sqr(vd) + 0.225*sqr(g1)*sqr(vu) + 0.375*sqr(g2)*sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);

   if (Mhh.minCoeff() < 0.)
      problems.flag_tachyon("hh");

   Mhh = sqrt(Mhh.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ah() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + sqr(Mu) + 0.3872983346207417*g1*g2*std::cos
      (ThetaW())*std::sin(ThetaW())*sqr(vd) + 0.075*sqr(g1)*sqr(vd) + 0.125*sqr(g2)*
      sqr(vd) - 0.075*sqr(g1)*sqr(vu) - 0.125*sqr(g2)*sqr(vu) + 0.25*sqr(g2)*
      sqr(vd)*sqr(std::cos(ThetaW())) + 0.15*sqr(g1)*sqr(vd)*sqr(std::sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.5*BMu + 0.5*BMu - 0.3872983346207417*g1*
      g2*vd*vu*std::cos(ThetaW())*std::sin(ThetaW()) - 0.25*vd*vu*sqr(g2)*sqr(std::cos(ThetaW(
      ))) - 0.15*vd*vu*sqr(g1)*sqr(std::sin(ThetaW()));
   mass_matrix_Ah(1,1) = mHu2 + sqr(Mu) - 0.075*sqr(g1)*sqr(vd) -
      0.125*sqr(g2)*sqr(vd) + 0.3872983346207417*g1*g2*std::cos(ThetaW())*std::sin(ThetaW
      ())*sqr(vu) + 0.075*sqr(g1)*sqr(vu) + 0.125*sqr(g2)*sqr(vu) + 0.25*sqr(g2
      )*sqr(vu)*sqr(std::cos(ThetaW())) + 0.15*sqr(g1)*sqr(vu)*sqr(std::sin(ThetaW()));

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);

   if (MAh.minCoeff() < 0.)
      problems.flag_tachyon("Ah");

   MAh = sqrt(MAh.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + sqr(Mu) + 0.075*sqr(g1)*sqr(vd) +
      0.375*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*sqr(g2)*sqr(vu);
   mass_matrix_Hpm(0,1) = BMu;
   mass_matrix_Hpm(1,1) = mHu2 + sqr(Mu) - 0.075*sqr(g1)*sqr(vd) +
      0.125*sqr(g2)*sqr(vd) + 0.075*sqr(g1)*sqr(vu) + 0.375*sqr(g2)*sqr(vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);

   if (MHpm.minCoeff() < 0.)
      problems.flag_tachyon("Hpm");

   MHpm = sqrt(MHpm.cwiseAbs());
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Chi() const
{
   Eigen::Matrix<double,4,4> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassB;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
   mass_matrix_Chi(2,2) = 0;
   mass_matrix_Chi(2,3) = -Mu;
   mass_matrix_Chi(3,3) = 0;

   Symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Cha;

   mass_matrix_Cha(0,0) = MassWB;
   mass_matrix_Cha(0,1) = 0.7071067811865475*g2*vu;
   mass_matrix_Cha(1,0) = 0.7071067811865475*g2*vd;
   mass_matrix_Cha(1,1) = Mu;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
}

double CLASSNAME::get_mass_matrix_VWm() const
{
   const double mass_matrix_VWm = 0.25*sqr(g2)*(sqr(vd) + sqr(vu));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{
   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = std::abs(mass_matrix_VWm);

   if (MVWm < 0.)
      problems.flag_tachyon("VWm");

   MVWm = std::sqrt(std::abs(MVWm));
}


double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = mHd2*vd + vd*sqr(Mu) - 0.5*vu*BMu - 0.5*vu*BMu +
      0.075*std::pow(vd,3)*sqr(g1) + 0.125*std::pow(vd,3)*sqr(g2) - 0.075*vd*sqr(g1)*
      sqr(vu) - 0.125*vd*sqr(g2)*sqr(vu);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   double result = mHu2*vu + vu*sqr(Mu) - 0.5*vd*BMu - 0.5*vd*BMu +
      0.075*std::pow(vu,3)*sqr(g1) + 0.125*std::pow(vu,3)*sqr(g2) - 0.075*vu*sqr(g1)*
      sqr(vd) - 0.125*vu*sqr(g2)*sqr(vd);

   return result;
}

double CLASSNAME::ThetaW() const
{
   return std::atan((0.7745966692414834*g1)/g2);
}

double CLASSNAME::v() const
{
   return 2*std::sqrt(sqr(MVWm)/sqr(g2));
}


std::ostream& operator<<(std::ostream& ostr, const MSSMNoFV_onshell_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace gm2calc
