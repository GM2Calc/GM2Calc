//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Loop functions
///  - flavour violating decays of charged leptons (from hep-ph/9403398)
///  - RK from 1706.07570
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2017 Aug, 2018 Feb
///
///  \author Cristian Sierra 
///          (cristian.sierra@monash.edu)
///  \date 2020 Oct
///
///  \author Douglas Jacob
///          (douglas.jacob@monash.edu)
///  \date 2020 Nov
///
///  *********************************************

#ifndef __loop_functions_hpp
#define __loop_functions_hpp

#include <complex>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_result.h>


namespace Old_Code
{

    void THDM_mumugamma(triplet<double> &result)

    //Yukawas vertices for THDM
    namespace Vertices_THDM
    {
      std::complex<double> yff_h(int f, int fp, double mf, Eigen::Matrix3cd xi_f, double vev, double cosab);
    
      std::complex<double> yff_H(int f, int fp, double mf, Eigen::Matrix3cd xi_f, double vev, double cosab);
    }

    //Auxiliary function to flip the sign and choose the scalar field in the fermion-fermion-Higgs vertex
    namespace Yukawas
    {
      std::complex<double> yff_phi(int f, int i, int j, int phi, double mf, Eigen::Matrix3cd xi_f, Eigen::Matrix3cd VCKM, double vev, double cosab);
    }

    //2-loop functions for BR(l>l' gamma) for the gTHDM
    namespace TwoLoopFunctions
    {
      std::complex<double> TwoLoopFH(std::complex<double> x);
      std::complex<double> TwoLoopFA(std::complex<double> x);
      std::complex<double> TwoLoopGW(std::complex<double> x);

      // Two Loop Functions for Muon g-2 for gTHDM
      std::complex<double> TwoLoopPhi(double m1, double m2, double m3);

      std::complex<double> TwoLoopfgammaphi(double Nc, double Qf, double alph, double mmu, double mf, double mphi, double mW, double mZ);
      std::complex<double> TwoLoopfZbosonphi(double Nc, double Qf, double alph, double mmu, double mf, double mphi, double mW, double mZ, double glv, double gfv);
      std::complex<double> TwoLoopfgammaA(double Nc, double Qf, double alph, double mmu, double mf, double mphi, double mW, double mZ);
      std::complex<double> TwoLoopfZbosonA(double Nc, double Qf, double alph, double mmu, double mf, double mphi, double mW, double mZ, double glv, double gfv);

      std::complex<double> TwoLoopfCl(int xl);
      std::complex<double> TwoLoopfCd(int xu, double xd, double Qu, double Qd);
      std::complex<double> TwoLoopfCu(int xu, double xd, double Qu, double Qd);

      std::complex<double> TwoLoopfC(int f, double Nc, double Qu, double Qd, double alph, double mmu, std::vector<double> mf, double mphi, double mW, double mZ);

      double TwoLoopf1(double x, void * params);
      double TwoLoopf2(double x, void * params);
      double TwoLoopf3(double x, void * params);
      double TwoLoopf4(double x, void * params);

      double TwoLoopF1(double w);
      double TwoLoopF2(double w);
      double TwoLoopF3(double w);
      double TwoLoopF4(double w);

      struct G_params {double wa; double wb; int n;};

      double TwoLoopg(double x, void * params);

      double TwoLoopG(double wa, double wb, int n);
    }

    // Loop functions for one loop diagrams
    namespace OneLoopFunctions
    {
      double OneLoopB(const double x);
      double OneLoopC(const double x);
      double OneLoopE(const double x);
      double OneLoopF(const double x);
    }

    // Amplitudes for BR(l>l' gamma) for the gTHDM
    namespace Amplitudes
    {
      //1-loop AL and AR amplitudes
      std::complex<double> A_loop1L(int f, int l, int li, int lp, int phi, std::vector<double> mnu, std::vector<double> ml, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double vev, double cosab);
      std::complex<double> A_loop1R(int f, int l, int li, int lp, int phi, std::vector<double> mnu, std::vector<double> ml, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double vev, double cosab);

      //2-loop fermionic contribution
      //AL
      std::complex<double> A_loop2fL(int fe, int lf, int l, int lp, int phi, double ml, double mlf, double mphi, double mZ, double Qf, double QfZ, Eigen::Matrix3cd xi_f, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab);
      //AR
      std::complex<double> A_loop2fR(int fe, int lf, int l, int lp, int phi, double ml, double mlf, double mphi, double mZ, double Qf, double QfZ, Eigen::Matrix3cd xi_f, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab);

      //2-loop bosonic contribution
      //AL
      std::complex<double> A_loop2bL(int f, int l, int lp, int phi, double ml, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab, double mW, double mZ);
      //AR
      std::complex<double> A_loop2bR(int f, int l, int lp, int phi, double ml, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab, double mW, double mZ);
    }

    // Two Loop Contributions for gTHDM from 1607.06292
    namespace TwoLoopContributions
    {
      std::complex<double> gm2mu_loop2f(int fe, int lf, int l, int lp, int phi, double mmu, std::vector<double> mf, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd xi_f, Eigen::Matrix3cd VCKM, int Nc, std::vector<double> Qf, std::vector<double> gfv, double vev, double cosab, double mW, double mZ, double alph);

      double gm2mu_barrzeephigammaf(int fe, int lf, int l, int lp, int phi, double mmu, double mf, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd xi_f, Eigen::Matrix3cd VCKM, int Nc, double Qf, double vev, double cosab, double alph);
      double gm2mu_barrzeephigammaC(int fe, int l, int lp, int phi, double mmu, double mHp, double mphi, double couplingphiCC, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double vev, double cosab, double alph);
      double gm2mu_barrzeephigammaW(int fe, int l, int lp,int phi, double mmu, double mW, double mphi, double couplingphiWW, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double vev, double cosab, double alph);

      std::complex<double> gm2mu_barrzeeCHiggsWBosontb(int fe, int l, int lp, double mmu, std::vector<double> mf, double mHp, std::vector<double> Qf, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd xi_D, Eigen::Matrix3cd xi_U, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab, double mW, double mZ, double alph);
      double gm2mu_barrzeeCHiggsWBosonC(int fe, int l, int lp, int phi, double mmu, double mphi, double mHp, double couplingphiCC, std::complex<double> couplingphiCW, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab, double mW, double mZ, double alph);
      double gm2mu_barrzeeCHiggsWBosonW(int fe, int l, int lp, int phi, double mmu, double mphi, double mHp, double couplingphiCC, std::complex<double> couplingphiCW, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab, double mW, double mZ, double alph);
    }
    

}

#endif //#defined __loop_functions_hpp__
