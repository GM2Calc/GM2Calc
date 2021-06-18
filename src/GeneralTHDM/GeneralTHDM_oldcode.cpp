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
///  \date 2020 Nov, Dec
///
///  *********************************************

namespace Old_Code
{

    // Calculation of muon g-2 
    void THDM_mumugamma(triplet<double> &result)
    {
      using namespace Pipes::THDM_mumugamma;

      #ifdef DEBUG
        cout<<"Starting THDM_mumugamma"<<endl;
      #endif

      SMInputs sminputs = *Dep::SMINPUTS;
      Spectrum spectrum = *Dep::THDM_spectrum;
      const int f = 0, l = 1, lp = 1;

      const double Alpha = 1/(sminputs.alphainv);
      const double alpha = spectrum.get(Par::dimensionless,"alpha");
      const double tanb = spectrum.get(Par::dimensionless,"tanb");
      const double beta = atan(tanb);
      const double cosb = cos(beta);
      const double v = sqrt(1.0/(sqrt(2.0)*sminputs.GF));
      const double cab = cos(alpha-beta);
      const double mW = sminputs.mW;
      const double mZ = sminputs.mZ;
      const double mE = sminputs.mE;
      const double mMu = sminputs.mMu;
      const double mTau = sminputs.mTau;
      const double mNu1 = sminputs.mNu1;
      const double mNu2 = sminputs.mNu2;
      const double mNu3 = sminputs.mNu3;
      const double mBmB = sminputs.mBmB;
      const double mT = sminputs.mT;
      const double mCmC = sminputs.mCmC;
      const double mS = sminputs.mS;
      const double m122 = spectrum.get(Par::mass1, "m12_2");
      const double mh = spectrum.get(Par::Pole_Mass,"h0",1);
      const double mH = spectrum.get(Par::Pole_Mass,"h0",2);
      const double mA = spectrum.get(Par::Pole_Mass,"A0");
      const double mHp = spectrum.get(Par::Pole_Mass,"H+");
      const std::vector<double> ml = {mE, mMu, mTau};     // charged leptons
      const std::vector<double> mvl = {mNu1, mNu2, mNu3}; // neutrinos
      const std::vector<double> mlf = {mTau, mBmB, mT};   // fermions in the second loop
      const std::vector<double> mphi = {mh, mH, mA, mHp};
      const double Yee = spectrum.get(Par::dimensionless,"Ye2",1,1);
      const double Yemu = spectrum.get(Par::dimensionless,"Ye2",1,2);
      const double Ymue = spectrum.get(Par::dimensionless,"Ye2",2,1);
      const double Yetau = spectrum.get(Par::dimensionless,"Ye2",1,3);
      const double Ytaue = spectrum.get(Par::dimensionless,"Ye2",3,1);
      const double Ymumu = spectrum.get(Par::dimensionless,"Ye2",2,2);
      const double Ymutau = spectrum.get(Par::dimensionless,"Ye2",2,3);
      const double Ytaumu = spectrum.get(Par::dimensionless,"Ye2",3,2);
      const double Ytautau = spectrum.get(Par::dimensionless,"Ye2",3,3);
      const double Ytt = spectrum.get(Par::dimensionless,"Yu2",3,3);
      const double Ytc = spectrum.get(Par::dimensionless,"Yu2",3,2);
      const double Ycc = spectrum.get(Par::dimensionless,"Yu2",2,2);
      const double Ybb = spectrum.get(Par::dimensionless,"Yd2",3,3);
      const double Ysb = spectrum.get(Par::dimensionless,"Yd2",2,3);
      const double Yss = spectrum.get(Par::dimensionless,"Yd2",2,2);
      const double A      = sminputs.CKM.A;
      const double lambda = sminputs.CKM.lambda;
      const double rhobar = sminputs.CKM.rhobar;
      const double etabar = sminputs.CKM.etabar;
      const std::complex<double> Vud(1 - (1/2)*lambda*lambda);
      const std::complex<double> Vcd(-lambda,0);
      const std::complex<double> Vtd((1-rhobar)*A*pow(lambda,3),-etabar*A*pow(lambda,3));
      const std::complex<double> Vus(lambda,0);
      const std::complex<double> Vcs(1 - (1/2)*lambda*lambda,0);
      const std::complex<double> Vts(-A*lambda*lambda,0);
      const std::complex<double> Vub(rhobar*A*pow(lambda,3),-etabar*A*pow(lambda,3));
      const std::complex<double> Vcb(A*lambda*lambda,0);
      const std::complex<double> Vtb(1,0);
      const double xitt = -((sqrt(2)*mT*tanb)/v) + Ytt/cosb;
      const double xicc = -((sqrt(2)*mCmC*tanb)/v) + Ycc/cosb;
      const double xitc = Ytc/cosb;
      const double xibb = -((sqrt(2)*mBmB*tanb)/v) + Ybb/cosb;
      const double xiss = -((sqrt(2)*mS*tanb)/v) + Yss/cosb;
      const double xisb = Ysb/cosb;
      const double xiee = -((sqrt(2)*mE*tanb)/v) + Yee/cosb;
      const double xiemu = Yemu/cosb;
      const double ximue = Ymue/cosb;
      const double xietau = Yetau/cosb;
      const double xitaue = Ytaue/cosb;
      const double ximumu = -((sqrt(2)*mMu*tanb)/v) + Ymumu/cosb;
      const double ximutau = Ymutau/cosb;
      const double xitaumu = Ytaumu/cosb;
      const double xitautau = -((sqrt(2)*mTau*tanb)/v) + Ytautau/cosb;

      Eigen::Matrix3cd xi_L, xi_U, xi_D, VCKM;

      xi_L << xiee,  xiemu,  xietau,
              ximue, ximumu, ximutau,
              xitaue, xitaumu, xitautau;

      xi_U << 0,   0,    0,
              0, xicc,  xitc,
              0, xitc, xitt;

      xi_D << 0,   0,    0,
              0, xiss,  xisb,
              0, xisb, xibb;

      const std::vector<Eigen::Matrix3cd> xi_f = {xi_L, xi_D, xi_U};

      // Needed for Hpm-l-vl couplings
      VCKM << Vud, Vus, Vub,
              Vcd, Vcs, Vcb,
              Vtd, Vts, Vtb;

      // One loop amplitude
      std::complex<double> Aloop1L = 0, A1L = 0;
      std::complex<double> Aloop1R = 0, A1R = 0;
      //Charged higgs contributions are being neglected
      //no longer
      for (int phi=0; phi<=3; ++phi)
      {
        for (int li = 0; li <=2; ++li)
        {
          //if ((phi == 0) and (li == 1))
          //{
          // Ignore purely SM contributions with SM higgs and no flavour changing 
          //  Aloop1L += 0.0;
          //  Aloop1R += 0.0;
          //}
          //else
          //{
            A1L=(ml[l]*ml[lp]/(16*pow(pi*mphi[phi],2)))*Amplitudes::A_loop1L(f, l, li, lp, phi, mvl, ml, mphi[phi], xi_L, VCKM, v, cab);
            A1R=(ml[l]*ml[lp]/(16*pow(pi*mphi[phi],2)))*Amplitudes::A_loop1R(f, l, li, lp, phi, mvl, ml, mphi[phi], xi_L, VCKM, v, cab);
            Aloop1L += A1L;
            Aloop1R += A1R;
          //}
        }
      }

      // Need to remove the SM Higgs contribution
      // Use lighter Higgs as SM Higgs, set cab=0 to simulate SM Yukawas
      // Alternatively could use: 1607.06292, eqn (32)
      std::complex<double> Aloop1SML = (ml[l]*ml[lp]/(16*pow(pi*mphi[0],2)))*Amplitudes::A_loop1L(f, l, l, lp, 0, mvl, ml, mphi[0], xi_L, VCKM, v, 0.);
      std::complex<double> Aloop1SMR = (ml[l]*ml[lp]/(16*pow(pi*mphi[0],2)))*Amplitudes::A_loop1L(f, l, l, lp, 0, mvl, ml, mphi[0], xi_L, VCKM, v, 0.);

      // Two loop amplitude
      const std::vector<double> Qf = {-1.,-1./3.,2./3.};
      const std::vector<double> Nc = { 1.,    3.,   3.};

      const double sw2 = 1 - pow(mW/mZ,2);
      const std::vector<double> gfv = {-1./2./2.-Qf[0]*sw2, -1./2./2.-Qf[1]*sw2, 1./2./2.-Qf[2]*sw2};

      //Fermionic contribution, source: 1607.06292
      std::complex<double> Aloop2f = 0.;
      std::complex<double> Aloop2SMf = 0.;
      for (int phi=0; phi<=3; ++phi)
      { 
        for (int lf=0; lf<=2; ++lf)
        {
          Aloop2f += TwoLoopContributions::gm2mu_loop2f(f, lf, l, lp, phi, mMu, mlf, mphi[phi], xi_L, xi_f[lf], VCKM, Nc[lf], Qf, gfv, v, cab, mW, mZ, Alpha);
        }
      }

      // Use lighter Higgs as SM Higgs, set cab=0 to simulate SM Yukawas
      for (int lf=0; lf<=2; ++lf)
      {
        Aloop2SMf += TwoLoopContributions::gm2mu_loop2f(f, lf, l, lp, 0, mMu, mlf, mphi[0], xi_L, xi_f[lf], VCKM, Nc[lf], Qf, gfv, v, 0., mW, mZ, Alpha);
      }

      const std::vector<double> couplingphiCC = { \
      (-(mh*mh-2.*mHp*mHp) * cos(alpha-3.*beta) * sin(2.*beta) + cos(alpha+beta) * (-8.*m122+(3.*pow(mh,2)+2.*pow(mHp,2))*sin(2.*beta))) / pow(cos(beta)*sin(beta),2) / (8.*v*v), \
      (-(mH*mH-2.*mHp*mHp) * sin(alpha-3.*beta) + (3.*pow(mH,2)+2.*pow(mHp,2)-4.*m122/sin(beta)/cos(beta)) * sin(alpha + beta)) / sin(2.*beta) / (2.*v*v), \
      0.};
      const std::vector<double> couplingphiWW = {sqrt(1-pow(cab,2)), cab, 0.};
      const std::vector<std::complex<double>> couplingphiCW = { std::complex<double> (cab,0.),  std::complex<double> (-sqrt(1-pow(cab,2)),0.), std::complex<double> (0.,-1.)};
      
      //Barr-Zee contribution, source: 1502.04199
      std::complex<double> Aloop2BZ = 0.;
      std::complex<double> Aloop2SMBZ = 0.;
      for (int phi=0; phi<=2; ++phi)
      { 
        // Superseded by gm2mu_loop2f by neutral boson contributions
        //for (int lf=0; lf<=2; ++lf)
        //{
        //  Aloop2BZ += TwoLoopContributions::gm2mu_barrzeephigammaf(f, lf, l, lp, phi, mMu, mlf[f], mphi[phi], xi_L, xi_f[lf], VCKM, Nc[f], Qf[f], v, cab, Alpha);
        //}
        Aloop2BZ += TwoLoopContributions::gm2mu_barrzeephigammaC(f, l, lp, phi, mMu, mphi[3], mphi[phi], couplingphiCC[phi], xi_L, VCKM, v, cab, Alpha);
        Aloop2BZ += TwoLoopContributions::gm2mu_barrzeephigammaW(f, l, lp, phi, mMu, mW, mphi[phi], couplingphiWW[phi], xi_L, VCKM, v, cab, Alpha);
        Aloop2BZ += TwoLoopContributions::gm2mu_barrzeeCHiggsWBosonC(f, l, lp, 3, mMu, mphi[3], mphi[phi], couplingphiCC[phi], couplingphiCW[phi], xi_L, VCKM, sw2, v, cab, mW, mZ, Alpha);
        Aloop2BZ += TwoLoopContributions::gm2mu_barrzeeCHiggsWBosonW(f, l, lp, 3, mMu, mphi[3], mphi[phi], couplingphiWW[phi], couplingphiCW[phi], xi_L, VCKM, sw2, v, cab, mW, mZ, Alpha);
      }
      // Superseded by gm2mu_loop2f by charged boson contributions
      //Aloop2BZ += TwoLoopContributions::gm2mu_barrzeeCHiggsWBosontb(f, l, lp, mMu, mlf, mHp, Qf, xi_L, xi_D, xi_U, VCKM, sw2, v, cab, mW, mZ, Alpha);

      // Use lighter Higgs as SM Higgs, set cab=0 to simulate SM Yukawas
      Aloop2SMBZ += TwoLoopContributions::gm2mu_barrzeephigammaW(f, l, lp, 0, mMu, mW, mphi[0], couplingphiWW[0], xi_L, VCKM, v, 0., Alpha);

      //Bosonic contribution
      // 3-boson contributions suppressed and neglected

      result.central = Aloop1L.real() + Aloop1R.real() - Aloop1SML.real() - Aloop1SMR.real() + Aloop2f.real() - Aloop2SMf.real() + Aloop2BZ.real() - Aloop2SMBZ.real();
      result.upper = std::max(std::abs(result.central)*0.3, 6e-10); //Based on hep-ph/0609168v1 eqs 84 & 85
      result.lower = result.upper;

      #ifdef DEBUG
        printf("(g-2)_mu=%.3e\n",result.central);
        cout<<"Finished THDM_mumugamma"<<endl;
      #endif
    }

  //Yukawas vertices for THDM
  namespace Vertices_THDM
  {
    std::complex<double> yff_h(int f, int fp, double mf, Eigen::Matrix3cd xi_f, double vev, double cosab)
    {
      Eigen::Matrix3i delta_ij;
      delta_ij << 1.0,  0.0,  0.0,
                  0.0,  1.0,  0.0,
                  0.0,  0.0,  1.0;

      return (mf/vev)*(std::sqrt(1-cosab*cosab))*delta_ij(f,fp) + (cosab/std::sqrt(2))*xi_f(f,fp);
    }

    std::complex<double> yff_H(int f, int fp, double mf, Eigen::Matrix3cd xi_f, double vev, double cosab)
    {
      Eigen::Matrix3i delta_ij;
      delta_ij << 1.0,  0.0,  0.0,
                  0.0,  1.0,  0.0,
                  0.0,  0.0,  1.0;

      return (mf/vev)*cosab*delta_ij(f,fp) - (std::sqrt(1-cosab*cosab)/std::sqrt(2))*xi_f(f,fp);
    }
  }

  //Auxiliary function to flip the sign and choose the scalar field in the fermion-fermion-Higgs vertex
  namespace Yukawas
  {
    std::complex<double> yff_phi(int f, int i, int j, int phi, double mf, Eigen::Matrix3cd xi_f, Eigen::Matrix3cd VCKM, double vev, double cosab)
    {
       std::complex<double> I(0,1);
       switch(phi)
       {
         case 0 :
           return Vertices_THDM::yff_h(i, j, mf, xi_f, vev, cosab);
         case 1 :
           return Vertices_THDM::yff_H(i, j, mf, xi_f, vev, cosab);
         case 2 :
           if (f==2)//flip sign if it is an up-like quark
           {
             return -I*(1/std::sqrt(2))*xi_f(i,j);
           }
           else
           {
             return I*(1/std::sqrt(2))*xi_f(i,j);
           }
         case 3 :
           std::complex<double> vertex_total(0,0);
           switch(f)
           {
             case 2 ://flip sign and include CKM Matrix if it is an up-like quark
               for (int kk = 0; kk <= 2; ++kk)
               {
                 vertex_total += -VCKM(i,kk)*xi_f(kk,j);
               }
               break;
             case 1 ://include CKM Matrix if it is an down-like quark
               for (int kk = 0; kk <= 2; ++kk)
               {
                 vertex_total += VCKM(i,kk)*xi_f(kk,j);
               }
               break;
             case 0 ://if it is a charged lepton
               vertex_total += xi_f(i,j);
               break;
           }
           return vertex_total;
       }
    }
  }

    //2-loop functions for BR(l>l' gamma) for the gTHDM
    namespace TwoLoopFunctions
    {
        std::complex<double> TwoLoopFH(std::complex<double> x)
        {
           std::complex<double> z(1 - 4*real(x),0);
           if (imag(std::sqrt(z))!=0 or std::real(std::sqrt(z))>0)
           {
             std::complex<double> one(1,0);
             std::complex<double> two(2,0);
             std::complex<double> z(1 - 4*real(x),0);
             std::complex<double> w1 = (-one + std::sqrt(z))/(one + std::sqrt(z));
             std::complex<double> w2 = (one + std::sqrt(z))/(-one + std::sqrt(z));
             std::complex<double> Logs = (-x/(two*std::sqrt(z)))*(-two*two*std::sqrt(z) - two*std::sqrt(z)*std::log(x) - std::log(w1)*std::log(x) + two*x*std::log(w1)*std::log(x) + std::log(w2)*std::log(x) - two*x*std::log(w2)*std::log(x));
             std::complex<double> w4 = (-two)/(-one + std::sqrt(z));
             std::complex<double> w5 = (two)/(one + std::sqrt(z));
             gsl_sf_result reD1, imD1;
             gsl_sf_result reD2, imD2;
             gsl_sf_complex_dilog_e(abs(w4), arg(w4), &reD1, &imD1);
             std::complex<double> Dilog1(reD1.val,imD1.val);
             gsl_sf_complex_dilog_e(abs(w5), arg(w5), &reD2, &imD2);
             std::complex<double> Dilog2(reD2.val,imD2.val);
             std::complex<double> Dilogs = (-x/(two*std::sqrt(z)))*((-two + two*two*x)*Dilog1 + (two - two*two*x)*Dilog2);
             // cout<<"OneLoopFH(x) = " << std::real(Logs+Dilogs) << endl;
             return std::real(Logs+Dilogs);
           }
           else
           {
             utils_error().raise(LOCAL_INFO, "1/0 in dilog OneLoopFH");
           }
         }

        std::complex<double> TwoLoopFA(std::complex<double> x)
        {
           std::complex<double> z(1 - 4*real(x),0);
           if (imag(std::sqrt(z))!=0 or std::real(std::sqrt(z))>0)
            {
             std::complex<double> one(1,0);
             std::complex<double> two(2,0);
             std::complex<double> z(1 - 4*real(x),0);
             std::complex<double> w1 = (-one + std::sqrt(z))/(one + std::sqrt(z));
             std::complex<double> w2 = (one + std::sqrt(z))/(-one + std::sqrt(z));
             std::complex<double> Logs = (-x/(two*std::sqrt(z)))*(-std::log(w1)+std::log(w2))*std::log(x);
             std::complex<double> w4 = (-two)/(-one + std::sqrt(z));
             std::complex<double> w5 = (two)/(one + std::sqrt(z));
             gsl_sf_result reD1, imD1;
             gsl_sf_result reD2, imD2;
             gsl_sf_complex_dilog_e(abs(w4), arg(w4), &reD1, &imD1);
             std::complex<double> Dilog1(reD1.val,imD1.val);
             gsl_sf_complex_dilog_e(abs(w5), arg(w5), &reD2, &imD2);
             std::complex<double> Dilog2(reD2.val,imD2.val);
             std::complex<double> Dilogs = (-x/(two*std::sqrt(z)))*(-two*Dilog1 + two*Dilog2);
             return std::real(Logs+Dilogs);
            }
            else
            {
              utils_error().raise(LOCAL_INFO, "1/0 in dilog TwoLoopFA");
            }
        }

        std::complex<double> TwoLoopGW(std::complex<double> x)
        {
           std::complex<double> z(1 - 4*real(x),0);
           if (imag(std::sqrt(z))!=0 or std::real(std::sqrt(z))>0)
           {
             std::complex<double> one(1,0);
             std::complex<double> two(2,0);
             std::complex<double> four(4,0);
             std::complex<double> z(1 - 4*real(x),0);
             std::complex<double> w1 = (-one + std::sqrt(z))/(one + std::sqrt(z));
             std::complex<double> w2 = (one + std::sqrt(z))/(-one + std::sqrt(z));
             std::complex<double> Logs =  - std::sqrt(-z)*std::log(-one - std::sqrt(z)) + four*x*std::sqrt(-z)*std::log(-one - std::sqrt(z))
                             + std::sqrt(-z)*std::log(one - std::sqrt(z)) - four*x*std::sqrt(-z)*std::log(one - std::sqrt(z))         
                             + std::sqrt(-z)*std::log(-one + std::sqrt(z)) - four*x*std::sqrt(-z)*std::log(-one + std::sqrt(z)) - std::sqrt(-z)*std::log(one + std::sqrt(z))
                             + four*x*std::sqrt(-z)*std::log(one + std::sqrt(z)) - two*std::sqrt(-std::pow(-z,2))*std::log(x)         
                             + two*x*std::sqrt(-z)*std::log(w1)*std::log(one/x) - two*x*std::sqrt(-z)*std::log(w2)*std::log(one/x);
             std::complex<double> w4 = (-two)/(-one + std::sqrt(z));
             std::complex<double> w5 = (two)/(one + std::sqrt(z));
             gsl_sf_result reD1, imD1;
             gsl_sf_result reD2, imD2;
             gsl_sf_complex_dilog_e(abs(w4), arg(w4), &reD1, &imD1);
             std::complex<double> Dilog1(reD1.val,imD1.val);
             gsl_sf_complex_dilog_e(abs(w5), arg(w5), &reD2, &imD2);
             std::complex<double> Dilog2(reD2.val,imD2.val);
             std::complex<double> Dilogs = four*x*std::sqrt(-z)*(-Dilog1 + Dilog2);
             std::complex<double> Atans = -four*std::sqrt(z)*atan(one/std::sqrt(-z))*(one-four*x);
             return std::real((x/(two*std::sqrt(z)*std::pow(-z,1.5)))*(Logs+Dilogs+Atans));
           }
           else
           {
             utils_error().raise(LOCAL_INFO, "1/0 in dilog TwoLoopGW");
           }
        }

      // Two Loop Functions for Muon g-2 for gTHDM
      //  Source:  1607.06292, eqn (68)
      std::complex<double> TwoLoopPhi(double m1, double m2, double m3)
      { 
        double mm1, mm2, mm3;
        // Order the masses in increasing order
        if (m1 <= m2)
        {
          if (m1 <= m3)
          {
            if (m2 <= m3)
            {
              mm1 = m1; mm2 = m2; mm3 = m3;
            }
            else
            {
              mm1 = m1; mm2 = m3; mm3 = m2;
            }
          }
          else
          {
            mm1 = m3; mm2 = m1; mm3 = m2;
          }
        }
        else
        {
          if (m1 <= m3)
          {
            mm1 = m2; mm2 = m1; mm3 = m3;
          }
          else
          {
            if (m2 <= m3)
            {              
              mm1 = m2; mm2 = m3; mm3 = m1;
            }
            else
            {
              mm1 = m3; mm2 = m2; mm3 = m1;
            }
          }
        }
        if (mm3 != 0)
        {
          std::complex<double> argum  = std::complex<double> (std::pow(mm1,4)+std::pow(mm2,4)+std::pow(mm3,4)-2*std::pow(mm1*mm2,2)-2*std::pow(mm2*mm3,2)-2*std::pow(mm3*mm1,2));
          std::complex<double> lambda = std::sqrt(argum);
          std::complex<double> alphap = (std::pow(mm3,2)+std::pow(mm1,2)-std::pow(mm2,2)-lambda) / (2*std::pow(mm3,2));
          std::complex<double> alpham = (std::pow(mm3,2)-std::pow(mm1,2)+std::pow(mm2,2)-lambda) / (2*std::pow(mm3,2));
          gsl_sf_result reD1, imD1;
          gsl_sf_complex_dilog_e(abs(alphap), arg(alphap), &reD1, &imD1);
          std::complex<double> Dilog1(reD1.val,imD1.val);
          gsl_sf_result reD2, imD2;
          gsl_sf_complex_dilog_e(abs(alpham), arg(alpham), &reD2, &imD2);
          std::complex<double> Dilog2(reD2.val,imD2.val);
          return lambda/2. * (2.*std::log(alphap)*std::log(alpham) - std::log(std::pow(mm1/mm3,2))*std::log(pow(mm2/mm3,2)) - 2.*Dilog1 - 2.*Dilog2 + std::pow(M_PI,2)/3.);
        } 
        else 
        { 
          utils_error().raise(LOCAL_INFO, "1/0 in dilog TwoLoopPhi");
        }
      }   
        
      //  Source:  1607.06292, eqns (54-57)
      std::complex<double> TwoLoopfgammaphi(int Nc, double Qf, double alph, double mmu, double mf, double mphi, double mW, double mZ)
      {
        double sw2 = 1-std::pow(mW/mZ,2);
        std::complex<double> Fphi = -2. + std::log(std::pow(mphi/mf,2)) - (std::pow(mphi,2)-2.*std::pow(mf,2))/std::pow(mphi,2) * TwoLoopFunctions::TwoLoopPhi(mphi,mf,mf)/(std::pow(mphi,2)-4.*std::pow(mf,2));
        return Nc * std::pow(Qf * alph * mmu,2) / pow(2.*M_PI*mW,2) / sw2 * std::pow(mf/mphi,2) * Fphi;
      }

      std::complex<double> TwoLoopfZbosonphi(int Nc, double Qf, double alph, double mmu, double mf, double mphi, double mW, double mZ, double glv, double gfv)
      {
        double sw2 = 1-std::pow(mW/mZ,2);
        double cw2 = 1-sw2;
        std::complex<double> Fphi = -2. + std::log(std::pow(mphi/mf,2)) - (std::pow(mphi,2)-2.*std::pow(mf,2))/std::pow(mphi,2) * TwoLoopFunctions::TwoLoopPhi(mphi,mf,mf)/(std::pow(mphi,2)-4.*std::pow(mf,2));
        std::complex<double> FZ   = -2. + std::log(std::pow(mZ/mf,2)) - (std::pow(mZ,2)-2.*std::pow(mf,2))/std::pow(mZ,2) * TwoLoopFunctions::TwoLoopPhi(mZ,mf,mf)/(std::pow(mZ,2)-4.*std::pow(mf,2));
        return -Nc * Qf * glv * gfv * std::pow(alph * mmu,2) / std::pow(2.*M_PI*mW*sw2,2) / cw2 * std::pow(mf,2)/(std::pow(mphi,2)-std::pow(mZ,2)) * (Fphi-FZ);
      }

      std::complex<double> TwoLoopfgammaA(int Nc, double Qf, double alph, double mmu, double mf, double mphi, double mW, double mZ)
      {
        double sw2 = 1-std::pow(mW/mZ,2);
        std::complex<double> Fphi = TwoLoopFunctions::TwoLoopPhi(mphi,mf,mf)/(std::pow(mphi,2)-4.*std::pow(mf,2));
        return Nc * std::pow(Qf * alph * mmu,2) / pow(2.*M_PI*mW,2) / sw2 * std::pow(mf/mphi,2) * Fphi;
      }

      std::complex<double> TwoLoopfZbosonA(int Nc, double Qf, double alph, double mmu, double mf, double mphi, double mW, double mZ, double glv, double gfv)
      {
        double sw2 = 1-std::pow(mW/mZ,2);
        double cw2 = 1-sw2;
        std::complex<double> Fphi = TwoLoopFunctions::TwoLoopPhi(mphi,mf,mf)/(std::pow(mphi,2)-4.*std::pow(mf,2));
        std::complex<double> FZ   = TwoLoopFunctions::TwoLoopPhi(mZ,mf,mf)/(std::pow(mZ,2)-4.*std::pow(mf,2));
        return -Nc * Qf * glv * gfv * std::pow(alph * mmu,2) / std::pow(2.*M_PI*mW*sw2,2) / cw2 * std::pow(mf,2)/(std::pow(mphi,2)-std::pow(mZ,2)) * (Fphi-FZ);
      }

      //  Source:  1607.06292, eqns (60-62)
      std::complex<double> TwoLoopfCl(double xl)
      {
        std::complex<double> z = 1.-1./xl;
        if (z != 0.)
        {
          gsl_sf_result reD, imD;
          gsl_sf_complex_dilog_e(abs(z), arg(z), &reD, &imD);
          std::complex<double> Dilog(reD.val,imD.val);
          return xl + xl * (xl-1.) * (Dilog-std::pow(M_PI,2)/6.) + (xl-0.5)*std::log(xl);
        }
        else
        {
          utils_error().raise(LOCAL_INFO, "1/0 in dilog TwoLoopfCl");
        }
      }

      std::complex<double> TwoLoopfCd(double xu, double xd, double Qu, double Qd)
      {
        double y  = std::pow(xu-xd,2) - 2.*(xu+xd) + 1.;
        double s  = (Qu + Qd) / 4.;
        double c  = std::pow(xu-xd,2) - Qu*xu + Qd*xd;
        double cb = (xu-Qu)*xu - (xd+Qd)*xd;
        std::complex<double> z = 1.-xd/xu;
        if (z != 0.)
        {
          gsl_sf_result reD, imD;
          gsl_sf_complex_dilog_e(abs(z), arg(z), &reD, &imD);
          std::complex<double> Dilog(reD.val,imD.val);
          return -(xu-xd) + (cb/y-c*(xu-xd)/y) * TwoLoopFunctions::TwoLoopPhi(std::sqrt(xd),std::sqrt(xu),1.) \
                 + c * (Dilog - 0.5*std::log(xu)*std::log(xd/xu) /** TwoLoopFunctions::TwoLoopPhi(std::sqrt(xd),std::sqrt(xu),1.))*/) \
                 + (s+xd) * std::log(xd) + (s-xu)*std::log(xu);
        }
        else
        {
          utils_error().raise(LOCAL_INFO, "1/0 in dilog TwoLoopfCd");
        }      
      }

      std::complex<double> TwoLoopfCu(double xu, double xd, double Qu, double Qd)
      {
        double y  = std::pow(xu-xd,2) - 2.*(xu+xd) + 1.;
        double s  = (Qu + 2. + Qd + 2.) / 4.;
        double c  = std::pow(xu-xd,2) - (Qu+2.)*xu + (Qd+2.)*xd;
        double cb = (xu-Qu-2.)*xu - (xd+Qd+2.)*xd;
        std::complex<double> z = 1.-xd/xu;
        if (z != 0.)
        {
          gsl_sf_result reD, imD;
          gsl_sf_complex_dilog_e(abs(z), arg(z), &reD, &imD);
          std::complex<double> Dilog(reD.val,imD.val);
          return -(xu-xd) + (cb/y-c*(xu-xd)/y) * TwoLoopFunctions::TwoLoopPhi(std::sqrt(xd),std::sqrt(xu),1.) \
                 + c * (Dilog - 0.5*std::log(xu)*std::log(xd/xu) /** TwoLoopFunctions::TwoLoopPhi(std::sqrt(xd),std::sqrt(xu),1.)*/) \
                 + (s+xd) * std::log(xd) + (s-xu)*std::log(xu) \
                 - 4./3. * (xu-xd-1.)/y * TwoLoopFunctions::TwoLoopPhi(std::sqrt(xd),std::sqrt(xu),1.) - 1./3.*(std::pow(std::log(xd),2)-std::pow(std::log(xu),2));
        }
        else
        {
          utils_error().raise(LOCAL_INFO, "1/0 in dilog TwoLoopfCu");
        }
      }

      std::complex<double> TwoLoopfC(int lf, int Nc, double Qu, double Qd, double alph, double mmu, std::vector<double> mf, double mphi, double mW, double mZ)
      {
        double sw2 = 1-std::pow(mW/mZ,2);
        if (lf == 0)
        {
          const double xlC = std::pow(mf[lf]/mphi,2);
          const double xlW = std::pow(mf[lf]/mW,2);
          return std::pow((alph*mmu)*mf[lf]/(M_PI*mW*sw2),2) * Nc / (std::pow(mphi,2)-std::pow(mW,2)) / 32. * (TwoLoopFunctions::TwoLoopfCl(xlC) - TwoLoopFunctions::TwoLoopfCl(xlW));
        } 
        else if (lf == 1)
        {
          const double xuC = std::pow(mf[2]/mphi,2);
          const double xdC = std::pow(mf[1]/mphi,2);
          const double xuW = std::pow(mf[2]/mW,2);
          const double xdW = std::pow(mf[1]/mW,2);
          return std::pow((alph*mmu)*mf[lf]/(M_PI*mW*sw2),2) * Nc / (std::pow(mphi,2)-std::pow(mW,2)) / 32. * (TwoLoopFunctions::TwoLoopfCd(xuC, xdC, Qu, Qd) - TwoLoopFunctions::TwoLoopfCd(xuW, xdW, Qu, Qd));
        }
        else if (lf == 2)
        {          
          const double xuC = std::pow(mf[2]/mphi,2);
          const double xdC = std::pow(mf[1]/mphi,2);
          const double xuW = std::pow(mf[2]/mW,2);
          const double xdW = std::pow(mf[1]/mW,2);
          return std::pow((alph*mmu)*mf[lf]/(M_PI*mW*sw2),2) * Nc / (std::pow(mphi,2)-std::pow(mW,2)) / 32. * (TwoLoopFunctions::TwoLoopfCu(xuC, xdC, Qu, Qd) - TwoLoopFunctions::TwoLoopfCu(xuW, xdW, Qu, Qd));
        }
      }

      //  Source:  1502.04199, eqn (25-29)
      double TwoLoopf1(double x, void * params)
      {
        double w = *(double *) params;
        return w/2. * (2.*x*(1.-x)-1.) / (w-x*(1.-x)) * std::log(w/(x*(1.-x)));
      }

      double TwoLoopf2(double x, void * params)
      {
        double w = *(double *) params;
        return 1./2. * x*(x-1.) / (w-x*(1.-x)) * std::log(w/(x*(1.-x)));
      }

      double TwoLoopf3(double x, void * params)
      {
        double w = *(double *) params;
        return 1./2. * (x*w*(3.*x*(4.*x-1.)+10.)-x*(1.-x)) / (w-x*(1.-x)) * std::log(w/(x*(1.-x)));
      }

      double TwoLoopf4(double x, void * params)
      {
        double w = *(double *) params;
        return w/2. * 1. / (w-x*(1.-x)) * std::log(w/(x*(1.-x)));
      }
      
      double TwoLoopF1(double w)
      {
        // Integral of w/2 * integral((2*x*(1-x)-1) / (w-x*(1-x)) * log(w/(x*(1-x))),{x,0,1})
        double result, error;
        const int alloc = 1000;
        gsl_integration_workspace * work = gsl_integration_workspace_alloc (alloc);
        gsl_function fun;
        fun.function = &TwoLoopf1;
        fun.params = &w;
        gsl_integration_qags (&fun, 0, 1, 0, 1e-7, alloc, work, &result, &error);
        gsl_integration_workspace_free(work);

        return result;
      }

      double TwoLoopF2(double w)
      {
        // Integral of 1/2 * integral(x*(x-1) / (w-x*(1-x)) * log(w/(x*(1-x))),{x,0,1})
        double result, error;
        const int alloc = 1000;
        gsl_integration_workspace * work = gsl_integration_workspace_alloc (alloc);
        gsl_function fun;
        fun.function = &TwoLoopf2;
        fun.params = &w;
        gsl_integration_qags (&fun, 0, 1, 0, 1e-7, alloc, work, &result, &error);
        gsl_integration_workspace_free(work);

        return result;
      }

      double TwoLoopF3(double w)
      {
        // Integral of 1/2 * integral((x*w*(3*x*(4*x-1)+10)-x*(1-x)) / (w-x*(1.-x)) * log(w/(x*(1-x))),{x,0,1})
        double result, error;
        const int alloc = 1000;
        gsl_integration_workspace * work = gsl_integration_workspace_alloc (alloc);
        gsl_function fun;
        fun.function = &TwoLoopf3;
        fun.params = &w;
        gsl_integration_qags (&fun, 0, 1, 0, 1e-7, alloc, work, &result, &error);
        gsl_integration_workspace_free(work);

        return result;
      }

      double TwoLoopF4(double w)
      {
        // Integral of w/2 * integral(1 / (w-x*(1-x)) * log(w/(x*(1-x))),{x,0,1})
        double result, error;
        const int alloc = 1000;
        gsl_integration_workspace * work = gsl_integration_workspace_alloc (alloc);
        gsl_function fun;
        fun.function = &TwoLoopf4;
        fun.params = &w;
        gsl_integration_qags (&fun, 0, 1, 0, 1e-7, alloc, work, &result, &error);
        gsl_integration_workspace_free(work);

        return result;
      }

      double TwoLoopg(double x, void * params)
      {
        struct G_params * p = (struct G_params *) params;
        double wa = (p->wa);
        double wb = (p->wb);
        int    n  = (p->n);
        return std::pow(x,n) * std::log((wa*x+wb*(1.-x)) / (x*(1.-x))) / (x*(1.-x)-wa*x-wb*(1.-x));
      }

      double TwoLoopG(double wa, double wb, int n)
      {
        // Integral of std::pow(x,n) * log((wa*x+wb*(1-x)) / (x*(1-x))) / (x*(1-x)-wa*x-wb*(1-x))
        double result, error;
        const int alloc = 1000;
        gsl_integration_workspace * work = gsl_integration_workspace_alloc (alloc);

        gsl_function fun;
        struct G_params params = {wa, wb, n};
        fun.function = &TwoLoopg;
        fun.params = &params;

        gsl_integration_qags (&fun, 0, 1, 0, 1e-7, alloc, work, &result, &error);

        gsl_integration_workspace_free(work);

        return result;
      }

    }

    // Loop functions for one loop diagrams
    // Source: FlexibleSUSY v2 manual 1710.03760, eqns (40-43)
    namespace OneLoopFunctions
    {
      double OneLoopB(const double x)
      {
        if(x == 0)
          return 2.;
        else if(x == 1)
          return 1.;
        else if(std::isinf(x))
          return 0.;
        else
          return 2.*(1. - 6.*x + 3.*x*x + 2.*x*x*x - 6.*x*x*std::log(x))/(std::pow(1.-x,4));
      }

      double OneLoopC(const double x)
      {
        if(x == 0)
          return 3.;
        else if(x == 1)
          return 1.;
        else if(std::isinf(x))
          return 0.;
        else
          return 3.*(1. - 1.*x*x + 2.*x*std::log(x))/(std::pow(1.-x,3));
      }

      double OneLoopE(const double x)
      {
        if(x == 0)
          return 4.;
        else if(x == 1)
          return 1.;
        else if(std::isinf(x))
          return 0.;
        else
          return 2.*(2. + 3.*x - 6.*x*x + 1.*x*x*x + 6.*x*std::log(x))/(std::pow(1.-x,4));
      }

      double OneLoopF(const double x)
      {
        if(x == 0)
          return 1./0.;
        else if(x == 1)
          return 1.;
        else if(std::isinf(x))
          return 0.;
        else
          return 3.*(-3. + 4.*x - 1.*x*x - 2.*std::log(x))/(2.*std::pow(1.-x,3));
      }
    }

    // Amplitudes for BR(l>l' gamma) for the gTHDM
    namespace Amplitudes
    {
      //1-loop AL and AR amplitudes
      std::complex<double> A_loop1L(int f, int l, int li, int lp, int phi, std::vector<double> mnu, std::vector<double> ml, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double vev, double cosab)
      {
        // f = 0,1,2 for electron,down,up families for external fermion
        // l,li,lp = 0,1,2 are the generation numbers of incoming,internal,outgoing lepton
        // phi = 0,1,2,3 for h,H,A,H+
        // mnu,ml is mass of neutrino,lepton in loop
        // mphi is array of Higgs boson masses
        //General contribution
        if ((phi == 0) or (phi == 1) or (phi == 2))
        {
          //FFS diagram
          double x = std::pow(ml[li]/mphi,2);
          std::complex<double> term1(ml[l]*ml[l]/std::pow(ml[l],2)  * OneLoopFunctions::OneLoopE(x)/24.,0.);
          std::complex<double> term2(ml[l]*ml[lp]/std::pow(ml[l],2) * OneLoopFunctions::OneLoopE(x)/24.,0.);
          std::complex<double> term3(ml[l]*ml[li]/std::pow(ml[l],2) * OneLoopFunctions::OneLoopF(x)/3., 0.);
          return term1*conj(Yukawas::yff_phi(f,li,lp,phi,ml[li],xi_L,VCKM,vev,cosab))*Yukawas::yff_phi(f,li,l,phi,ml[li],xi_L,VCKM,vev,cosab) \
                 + term2*conj(Yukawas::yff_phi(f,l,li,phi,ml[l],xi_L,VCKM,vev,cosab)) *Yukawas::yff_phi(f,lp,li,phi,ml[lp],xi_L,VCKM,vev,cosab) \
                 + term3*conj(Yukawas::yff_phi(f,li,lp,phi,ml[li],xi_L,VCKM,vev,cosab))*conj(Yukawas::yff_phi(f,l,li,phi,ml[l],xi_L,VCKM,vev,cosab));
        }
        else if (phi == 3)
        {
          //SSF diagram
          double x = std::pow(mnu[li]/mphi,2);
          std::complex<double> term1(ml[l]/ml[l] * OneLoopFunctions::OneLoopB(x)/24.,0.);
          return term1*conj(Yukawas::yff_phi(f,li,lp,phi,ml[li],xi_L,VCKM,vev,cosab))*Yukawas::yff_phi(f,li,l,phi,ml[li],xi_L,VCKM,vev,cosab);
        }
      }

      std::complex<double> A_loop1R(int f, int l, int li, int lp, int phi, std::vector<double> mnu, std::vector<double> ml, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double vev, double cosab)
      {
        // f = 0,1,2 for electron,down,up families for external fermion
        // l,li,lp = 0,1,2 are the generation numbers of incoming,internal,outgoing lepton
        // phi = 0,1,2,3 for h,H,A,H+
        // mnu,ml is mass of neutrino,lepton in loop
        // mphi is array of Higgs boson masses
        //General contribution
        if ((phi == 0) or (phi == 1) or (phi == 2))
        {
          //FFS diagram
          double x = std::pow(ml[li]/mphi,2);
          std::complex<double> term1(ml[l]*ml[l]/std::pow(ml[l],2)  * OneLoopFunctions::OneLoopE(x)/24.,0.);
          std::complex<double> term2(ml[l]*ml[lp]/std::pow(ml[l],2) * OneLoopFunctions::OneLoopE(x)/24.,0.);
          std::complex<double> term3(ml[l]*ml[li]/std::pow(ml[l],2) * OneLoopFunctions::OneLoopF(x)/3., 0.);
          return term1*conj(Yukawas::yff_phi(f,l,li,phi,ml[l],xi_L,VCKM,vev,cosab))*Yukawas::yff_phi(f,lp,li,phi,ml[lp],xi_L,VCKM,vev,cosab) \
                 + term2*conj(Yukawas::yff_phi(f,li,lp,phi,ml[li],xi_L,VCKM,vev,cosab)) *Yukawas::yff_phi(f,li,l,phi,ml[li],xi_L,VCKM,vev,cosab) \
                 + term3*Yukawas::yff_phi(f,li,l,phi,ml[li],xi_L,VCKM,vev,cosab)*Yukawas::yff_phi(f,lp,li,phi,ml[lp],xi_L,VCKM,vev,cosab);
        }
        else if (phi == 3)
        {
          //SSF diagram
          double x = std::pow(mnu[l]/mphi,2);
          std::complex<double> term1(ml[lp]/ml[l] * OneLoopFunctions::OneLoopB(x)/24.,0.);
          return term1*conj(Yukawas::yff_phi(f,li,lp,phi,ml[li],xi_L,VCKM,vev,cosab))*Yukawas::yff_phi(f,li,l,phi,ml[li],xi_L,VCKM,vev,cosab);
        }
      }

    //2-loop fermionic contribution
    //AL
    std::complex<double> A_loop2fL(int fe, int lf, int l, int lp, int phi, double ml, double mlf, double mphi, double mZ, double Qf, double QfZ, Eigen::Matrix3cd xi_f, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab)
    {
        std::complex<double> I(0,1);
        const int gi = 2; // Internal generation is the heaviest
        double xfphi = std::pow(mlf/mphi,2);
        double xfZ   = std::pow(mlf/mZ,2);
        std::complex<double> FHm = (xfphi*TwoLoopFunctions::TwoLoopFH(xfZ)-xfZ*TwoLoopFunctions::TwoLoopFH(xfphi))/(xfphi-xfZ);
        std::complex<double> FAm = (xfphi*TwoLoopFunctions::TwoLoopFA(xfZ)-xfZ*TwoLoopFunctions::TwoLoopFA(xfphi))/(xfphi-xfZ);
        return conj(Yukawas::yff_phi(fe, l, lp, phi, ml, xi_L, VCKM, vev, cosab))*(real(Yukawas::yff_phi(lf, gi, gi, phi, mlf, xi_f, VCKM, vev, cosab)) * (Qf*TwoLoopFunctions::TwoLoopFH(std::pow(mlf/mphi,2))+(1.-4.*sw2)*QfZ/(16.*sw2*(1-sw2))*FHm)-I*imag(Yukawas::yff_phi(lf, gi, gi, phi, mlf, xi_f, VCKM, vev, cosab)) * (Qf*TwoLoopFunctions::TwoLoopFA(std::pow(mlf/mphi,2))+(1.-4.*sw2)*QfZ/(16.*sw2*(1-sw2))*FAm));
    }
    
    //AR
    std::complex<double> A_loop2fR(int fe, int lf, int l, int lp, int phi, double ml, double mlf, double mphi, double mZ, double Qf, double QfZ, Eigen::Matrix3cd xi_f, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab)
    {
        std::complex<double> I(0,1);
        const int gi = 2; // Internal generation is the heaviest
        double xfphi = std::pow(mlf/mphi,2);
        double xfZ   = std::pow(mlf/mZ,2);
        std::complex<double> FHm = (xfphi*TwoLoopFunctions::TwoLoopFH(xfZ)-xfZ*TwoLoopFunctions::TwoLoopFH(xfphi))/(xfphi-xfZ);
        std::complex<double> FAm = (xfphi*TwoLoopFunctions::TwoLoopFA(xfZ)-xfZ*TwoLoopFunctions::TwoLoopFA(xfphi))/(xfphi-xfZ);
        return Yukawas::yff_phi(fe, l, lp, phi, ml, xi_L, VCKM, vev, cosab)*(real(Yukawas::yff_phi(lf, gi, gi, phi, mlf, xi_f, VCKM, vev, cosab)) * (Qf*TwoLoopFunctions::TwoLoopFH(std::pow(mlf/mphi,2))+(1.-4.*sw2)*QfZ/(16.*sw2*(1-sw2))*FHm)-I*imag(Yukawas::yff_phi(lf, gi, gi, phi, mlf, xi_f, VCKM, vev, cosab)) * (Qf*TwoLoopFunctions::TwoLoopFA(std::pow(mlf/mphi,2))+(1.-4.*sw2)*QfZ/(16.*sw2*(1-sw2))*FAm));
    }

      //2-loop bosonic contribution
      //AL
      std::complex<double> A_loop2bL(int f, int l, int lp, int phi, double ml, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab, double mW, double mZ)
      {
        double tw2 = sw2/(1-sw2);
        std::complex<double> xWphi(std::pow(mW/mphi,2),0);
        std::complex<double> xWZ(std::pow(mW/mZ,2),0);
        std::complex<double> FHm = (xWphi*TwoLoopFunctions::TwoLoopFH(xWZ)-xWZ*TwoLoopFunctions::TwoLoopFH(xWphi))/(xWphi-xWZ);
        std::complex<double> FAm = (xWphi*TwoLoopFunctions::TwoLoopFA(xWZ)-xWZ*TwoLoopFunctions::TwoLoopFA(xWphi))/(xWphi-xWZ);
        std::complex<double> pH(5-tw2+(1-tw2)/(2*std::pow(mW/mphi,2)),0);
        std::complex<double> pA(7-3*tw2-(1-tw2)/(2*std::pow(mW/mphi,2)),0);
        std::complex<double> onehalf(0.5,0);
        std::complex<double> two(2,0);
        std::complex<double> three(3,0);
        std::complex<double> four(4,0);
        std::complex<double> twenty3(23,0);
        return  conj(Yukawas::yff_phi(f, l, lp, phi, ml, xi_L, VCKM, vev, cosab))*(three*TwoLoopFunctions::TwoLoopFH(xWphi)+(twenty3/four)*TwoLoopFunctions::TwoLoopFA(xWphi)+(three/four)*TwoLoopFunctions::TwoLoopGW(xWphi)+(onehalf)*std::pow(mphi/mW,2)*(TwoLoopFunctions::TwoLoopFH(xWphi)-TwoLoopFunctions::TwoLoopFA(xWphi))+((1-4*sw2)/(8*sw2))*(pH*FHm + pA*FAm + (three/two)*(TwoLoopFunctions::TwoLoopFA(xWphi)+TwoLoopFunctions::TwoLoopGW(xWphi))));
      }
      //AR
      std::complex<double> A_loop2bR(int f, int l, int lp, int phi, double ml, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab, double mW, double mZ)
      {
        double tw2 = sw2/(1-sw2);
        std::complex<double> xWphi(std::pow(mW/mphi,2),0);
        std::complex<double> xWZ(std::pow(mW/mZ,2),0);
        std::complex<double> FHm = (xWphi*TwoLoopFunctions::TwoLoopFH(xWZ)-xWZ*TwoLoopFunctions::TwoLoopFH(xWphi))/(xWphi-xWZ);
        std::complex<double> FAm = (xWphi*TwoLoopFunctions::TwoLoopFA(xWZ)-xWZ*TwoLoopFunctions::TwoLoopFA(xWphi))/(xWphi-xWZ);
        std::complex<double> pH(5-tw2+(1-tw2)/(2*std::pow(mW/mphi,2)),0);
        std::complex<double> pA(7-3*tw2-(1-tw2)/(2*std::pow(mW/mphi,2)),0);
        std::complex<double> onehalf(0.5,0);
        std::complex<double> two(2,0);
        std::complex<double> three(3,0);
        std::complex<double> four(4,0);
        std::complex<double> twenty3(23,0);
        return  Yukawas::yff_phi(f, l, lp, phi, ml, xi_L, VCKM, vev, cosab)*(three*TwoLoopFunctions::TwoLoopFH(xWphi)+(twenty3/four)*TwoLoopFunctions::TwoLoopFA(xWphi)+(three/four)*TwoLoopFunctions::TwoLoopGW(xWphi)+(onehalf)*std::pow(mphi/mW,2)*(TwoLoopFunctions::TwoLoopFH(xWphi)-TwoLoopFunctions::TwoLoopFA(xWphi))+((1-4*sw2)/(8*sw2))*(pH*FHm + pA*FAm + (three/two)*(TwoLoopFunctions::TwoLoopFA(xWphi)+TwoLoopFunctions::TwoLoopGW(xWphi))));
      }
    }

    // Two Loop muon g-2 Contributions for gTHDM
    namespace TwoLoopContributions
    {
      // Source: 1607.06292, eqns (53,58)
      std::complex<double> gm2mu_loop2f(int fe, int lf, int l, int lp, int phi, double mmu, std::vector<double> mf, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd xi_f, Eigen::Matrix3cd VCKM, int Nc, std::vector<double> Qf, std::vector<double> gfv, double vev, double cosab, double mW, double mZ, double alph)
      { 
        // fe = 0,1,2 for electron,down,up families for external fermion
        // l,lf,lp = 0,1,2 are the generation numbers of incoming,loop,outgoing lepton
        // phi = 0,1,2,3 for h,H,A,H+
        // mf is array of loop fermion masses
        // mphi is array of Higgs boson masses
        // Qf is array of fermion electric charges
        // gfv is array of fermion Z-charges
        const int gi = 2; // Internal generation is the heaviest
        if ((phi == 0) or (phi == 1))
        {
          return (TwoLoopFunctions::TwoLoopfgammaphi(Nc, Qf[lf], alph, mmu, mf[lf], mphi, mW, mZ) + TwoLoopFunctions::TwoLoopfZbosonphi(Nc, Qf[lf], alph, mmu, mf[lf], mphi, mW, mZ, gfv[0], gfv[lf])) * vev / mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab) * vev / mf[lf] * Yukawas::yff_phi(lf, gi, gi, phi, mf[lf], xi_f, VCKM, vev, cosab);
        }
        else if (phi == 2)
        {
          return (TwoLoopFunctions::TwoLoopfgammaA(Nc, Qf[lf], alph, mmu, mf[lf], mphi, mW, mZ) + TwoLoopFunctions::TwoLoopfZbosonA(Nc, Qf[lf], alph, mmu, mf[lf], mphi, mW, mZ, gfv[0], gfv[lf])) * vev / mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab) * vev / mf[lf] * Yukawas::yff_phi(lf, gi, gi, phi, mf[lf], xi_f, VCKM, vev, cosab);
        }
        else if (phi == 3)
        {
          const double Qu = Qf[2], Qd = Qf[1];
          const std::complex<double> I(0,1);
          return TwoLoopFunctions::TwoLoopfC(lf, Nc, Qu, Qd, alph, mmu, mf, mphi, mW, mZ) * I * vev / mmu * Yukawas::yff_phi(fe, l, lp, 2, mmu, xi_L, VCKM, vev, cosab) * I * vev / mf[lf] * Yukawas::yff_phi(lf, gi, gi, 2, mf[lf], xi_f, VCKM, vev, cosab);
        }
      }

      // Source: 1502.04199, eqns (19-24)
      double gm2mu_barrzeephigammaf(int fe, int lf, int l, int lp, int phi, double mmu, double mf, double mphi, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd xi_f, Eigen::Matrix3cd VCKM, int Nc, double Qf, double vev, double cosab, double alph)
      {
        // fe = 0,1,2 for electron,down,up families for external fermion
        // l,lf,lp = 0,1,2 are the generation numbers of incoming,loop,outgoing lepton
        // phi = 0,1,2,3 for h,H,A,H+
        // mf is loop fermion mass
        // mphi is neutral scalar boson mass
        // Qf is loop fermion charge
        // Nc is loop fermion colour number
        const double x = std::pow(mf/mphi,2);
        const int gi = 2; // Internal generation is the heaviest
        double term1 = vev / mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab).real() * vev / mf * Yukawas::yff_phi(lf, gi, gi, phi, mf, xi_f, VCKM, vev, cosab).real() * TwoLoopFunctions::TwoLoopF1(x);
        double term2 = vev / mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab).imag() * vev / mf * Yukawas::yff_phi(lf, gi, gi, phi, mf, xi_f, VCKM, vev, cosab).imag() * TwoLoopFunctions::TwoLoopF4(x);
        return alph * Nc * std::pow(mmu * Qf / vev,2) / (4. * std::pow(M_PI,3)) * (term1 + term2);
      }

      double gm2mu_barrzeephigammaC(int fe, int l, int lp, int phi, double mmu, double mHp, double mphi, double couplingphiCC, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double vev, double cosab, double alph)
      {
        // fe = 0,1,2 for electron,down,up families for external fermion
        // l,lp = 0,1,2 are the generation numbers of incoming,outgoing lepton
        // phi = 0,1,2,3 for h,H,A,H+
        // mphi is neutral scalar boson mass
        // couplingphiCC is the coupling between scalar phi and two H+
        const double x = std::pow(mHp/mphi,2);
        return alph * std::pow(mmu/mphi,2) / (8. * std::pow(M_PI,3)) * vev / mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab).real() * couplingphiCC * TwoLoopFunctions::TwoLoopF2(x);
      }
        
      double gm2mu_barrzeephigammaW(int fe, int l, int lp, int phi, double mmu, double mW, double mphi, double couplingphiWW, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double vev, double cosab, double alph)
      {
        // fe = 0,1,2 for electron,down,up families for external fermion
        // l,lp = 0,1,2 are the generation numbers of incoming,outgoing lepton
        // phi = 0,1,2,3 for h,H,A,H+
        // mphi is neutral scalar boson mass
        // couplingphiWW is the coupling between scalar phi and two W+
        const double x = std::pow(mW/mphi,2);
        return alph * std::pow(mmu/vev,2) / (8. * std::pow(M_PI,3)) * vev / mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab).real() * couplingphiWW * TwoLoopFunctions::TwoLoopF3(x);
      }

      std::complex<double> gm2mu_barrzeeCHiggsWBosontb(int fe, int l, int lp, double mmu, std::vector<double> mf, double mHp, std::vector<double> Qf, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd xi_D, Eigen::Matrix3cd xi_U, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab, double mW, double mZ, double alph)
      {
        // fe = 0,1,2 for electron,down,up families for external fermion
        // l,lp = 0,1,2 are the generation numbers of incoming,outgoing lepton
        // mf is the array of loop fermion masses
        // Qf is the array of loop fermion charges
        const int phi = 2; // Use CP-odd Higgs boson
        const int fu = 2, fd = 1; // External fermions are both leptons
        const int gi = 2; // Internal generation is the heaviest
        const int Nc = 3; // Since internal fermions are top+bottom, number of colours is 3
        const double mt = mf[2], mb = mf[1];
        const double xtC = std::pow(mt/mHp,2);
        const double xbC = std::pow(mb/mHp,2);
        const double xtW = std::pow(mt/mW,2);
        const double xbW = std::pow(mb/mW,2);
        const double Qu = Qf[2], Qd = Qf[1];
        // Top contributions
        double term1 = std::pow(mt,2)*Qu * (std::conj(vev/mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab)) * vev/mt * Yukawas::yff_phi(fu, gi, gi, phi, mt, xi_U, VCKM, vev, cosab)).real() * (TwoLoopFunctions::TwoLoopG(xtC,xbC,2)+TwoLoopFunctions::TwoLoopG(xtC,xbC,3)-TwoLoopFunctions::TwoLoopG(xtW,xbW,2)-TwoLoopFunctions::TwoLoopG(xtW,xbW,3));
        double term2 = std::pow(mt,2)*Qd * (std::conj(vev/mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab)) * vev/mt * Yukawas::yff_phi(fu, gi, gi, phi, mt, xi_U, VCKM, vev, cosab)).real() * (TwoLoopFunctions::TwoLoopG(xtC,xbC,1)-TwoLoopFunctions::TwoLoopG(xtC,xbC,3)-TwoLoopFunctions::TwoLoopG(xtW,xbW,1)+TwoLoopFunctions::TwoLoopG(xtW,xbW,3));
        // Bottom contributions
        double term3 = std::pow(mb,2)*Qu * (std::conj(vev/mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab)) * vev/mb * Yukawas::yff_phi(fd, gi, gi, phi, mb, xi_D, VCKM, vev, cosab)).real() * (TwoLoopFunctions::TwoLoopG(xtC,xbC,2)-TwoLoopFunctions::TwoLoopG(xtC,xbC,3)-TwoLoopFunctions::TwoLoopG(xtW,xbW,2)+TwoLoopFunctions::TwoLoopG(xtW,xbW,3));
        double term4 = std::pow(mb,2)*Qd * (std::conj(vev/mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab)) * vev/mb * Yukawas::yff_phi(fd, gi, gi, phi, mb, xi_D, VCKM, vev, cosab)).real() * (TwoLoopFunctions::TwoLoopG(xtC,xbC,1)-2.*TwoLoopFunctions::TwoLoopG(xtC,xbC,2)+TwoLoopFunctions::TwoLoopG(xtC,xbC,3)-TwoLoopFunctions::TwoLoopG(xtW,xbW,1)+2.*TwoLoopFunctions::TwoLoopG(xtW,xbW,2)-TwoLoopFunctions::TwoLoopG(xtW,xbW,3));
        return alph*Nc*norm(VCKM(2,1))*std::pow(mmu/vev,2) / (32.*std::pow(M_PI,3)*sw2) / (std::pow(mHp,2)-std::pow(mW,2)) * (term1 + term2 + term3 + term4);
      }

      double gm2mu_barrzeeCHiggsWBosonC(int fe, int l, int lp, int phi, double mmu, double mHp, double mphi, double couplingphiCC, std::complex<double> couplingphiCW, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab, double mW, double mZ, double alph)
      {
        // fe = 0,1,2 for electron,down,up families for external fermion
        // l,lp = 0,1,2 are the generation numbers of incoming,outgoing lepton
        // phi = 0,1,2,3 for h,H,A,H+
        // mphi is neutral scalar boson mass
        // couplingphiCC is the coupling between scalar phi and two H+
        // couplingphiCW is the coupling between scalar phi, H+, and W+
        const double xSC = std::pow(mphi/mHp, 2);
        const double xSW = std::pow(mphi/mW,  2);
        const double xCW = std::pow(mHp/mW,   2);
        double term1 = (TwoLoopFunctions::TwoLoopG(1.,xSC,3)-TwoLoopFunctions::TwoLoopG(xCW,xSW,3));
        double term2 = (TwoLoopFunctions::TwoLoopG(1.,xSC,2)-TwoLoopFunctions::TwoLoopG(xCW,xSW,2));
        return alph*std::pow(mmu,2) / (64.*std::pow(M_PI,3)*sw2) / (std::pow(mHp,2)-std::pow(mW,2)) * std::real(std::conj(vev / mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab)) * couplingphiCC * couplingphiCW) * (term1-term2);
      }

      double gm2mu_barrzeeCHiggsWBosonW(int fe, int l, int lp, int phi, double mmu, double mHp, double mphi, double couplingphiWW, std::complex<double> couplingphiCW, Eigen::Matrix3cd xi_L, Eigen::Matrix3cd VCKM, double sw2, double vev, double cosab, double mW, double mZ, double alph)
      {
        // fe = 0,1,2 for electron,down,up families for external fermion
        // l,lp = 0,1,2 are the generation numbers of incoming,outgoing lepton
        // phi = 0,1,2,3 for h,H,A,H+
        // mphi is neutral scalar boson mass
        // couplingphiWW is the coupling between scalar phi and two W+
        // couplingphiCW is the coupling between scalar phi, H+, and W+
        const double xSC = std::pow(mphi/mHp,2);
        const double xSW = std::pow(mphi/mW, 2);
        const double xWC = std::pow(mW/mHp,  2);
        double term1 = (std::pow(mHp,2)-3.*std::pow(mW,2)-std::pow(mphi,2)) * (TwoLoopFunctions::TwoLoopG(xWC,xSC,2)-TwoLoopFunctions::TwoLoopG(1.,xSW,2));
        double term2 = (std::pow(mHp,2)+   std::pow(mW,2)-std::pow(mphi,2)) * (TwoLoopFunctions::TwoLoopG(xWC,xSC,3)-TwoLoopFunctions::TwoLoopG(1.,xSW,3));
        return alph*std::pow(mmu/vev,2) / (64.*std::pow(M_PI,3)*sw2) / (std::pow(mHp,2)-std::pow(mW,2)) * std::real(std::conj(vev / mmu * Yukawas::yff_phi(fe, l, lp, phi, mmu, xi_L, VCKM, vev, cosab)) * couplingphiWW * couplingphiCW) * (term1-term2);
      }

    }


}


