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

#include "MSSMNoFV_onshell.hpp"
#include "numerics2.hpp"
#include "gm2_error.hpp"
#include "gm2_1loop.hpp"
#include "gm2_mb.hpp"
#include "ffunctions.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <algorithm>
#include <limits>

#include <boost/math/tools/roots.hpp>

#define WARNING(message)                                                \
   do { std::cerr << "Warning: " << message << '\n'; } while (0)

namespace gm2calc {

namespace {
   const double ALPHA_EM_THOMPSON = 1./137.035999074;
   const double DELTA_ALPHA_EM_MZ =
      + 0.031498 /*leptonic*/
      - 0.0000728 /*top*/
      + 0.027626 /*hadronic, arXiv:1105.3149v2 */;
   const double ALPHA_EM_MZ =
      ALPHA_EM_THOMPSON / (1. - DELTA_ALPHA_EM_MZ);

   double calculate_e(double alpha) {
      return std::sqrt(4. * Pi * alpha);
   }
   double calculate_alpha(double e) {
      return e * e / (4. * Pi);
   }
   /// returns a/b if result is finite, otherwise 0.
   double divide_finite(double a, double b) noexcept {
      const double result = a / b;
      return std::isfinite(result) ? result : 0.;
   }
   /// element-wise division
   template <int Rows, int Cols>
   Eigen::Matrix<double,Rows,Cols> cwise_div(
         const Eigen::Matrix<double,Rows,Cols>& a,
         const Eigen::Matrix<double,Rows,Cols>& b) noexcept
   {
      Eigen::Matrix<double,Rows,Cols> result(Eigen::Matrix<double,Rows,Cols>::Zero());

      for (int i = 0; i < Rows; i++) {
         for (int k = 0; k < Cols; k++) {
            result(i,k) = divide_finite(a(i,k), b(i,k));
         }
      }

      return result;
   }
} // anonymous namespace

MSSMNoFV_onshell::MSSMNoFV_onshell()
   : MSSMNoFV_onshell_mass_eigenstates()
   , verbose_output(false)
   , EL(calculate_e(ALPHA_EM_MZ))
   , EL0(calculate_e(ALPHA_EM_THOMPSON))
   , Ae(Eigen::Matrix<double,3,3>::Zero())
   , Au(Eigen::Matrix<double,3,3>::Zero())
   , Ad(Eigen::Matrix<double,3,3>::Zero())
{
   set_g3(calculate_e(0.1184));
   get_physical().MFt = 173.34;
   get_physical().MFb = 4.18;
   get_physical().MFe = 0.000510998928;
   get_physical().MFm  = 0.1056583715;
   get_physical().MFtau = 1.777;
   get_physical().MVWm = 80.385;
   get_physical().MVZ  = 91.1876;
   set_scale(get_physical().MVZ);
}

MSSMNoFV_onshell::MSSMNoFV_onshell(const MSSMNoFV_onshell_mass_eigenstates& model_)
   : MSSMNoFV_onshell_mass_eigenstates(model_)
   , verbose_output(false)
   , EL(calculate_e(ALPHA_EM_MZ))
   , EL0(calculate_e(ALPHA_EM_THOMPSON))
   , Ae(cwise_div(model_.get_TYe(), model_.get_Ye()))
   , Au(cwise_div(model_.get_TYu(), model_.get_Yu()))
   , Ad(cwise_div(model_.get_TYd(), model_.get_Yd()))
{
}

void MSSMNoFV_onshell::set_alpha_MZ(double alpha)
{
   EL = calculate_e(alpha);
}

void MSSMNoFV_onshell::set_alpha_thompson(double alpha)
{
   EL0 = calculate_e(alpha);
}

void MSSMNoFV_onshell::set_TB(double tanb)
{
   const double vev = get_vev();
   const double sinb = tanb / std::sqrt(1 + tanb*tanb);
   const double cosb = 1.   / std::sqrt(1 + tanb*tanb);
   set_vd(vev * cosb);
   set_vu(vev * sinb);
}

double MSSMNoFV_onshell::get_vev() const
{
   const double cW = get_MW() / get_MZ();
   const double g2 = get_EL() / std::sqrt(1. - cW*cW);
   const double vev = 2. * get_MW() / g2;

   return vev;
}

double MSSMNoFV_onshell::get_TB() const
{
   if (MSSMNoFV_onshell::is_zero(get_vd()))
      return 0.;
   return get_vu() / get_vd();
}

/**
 * Converts the model parameters from the DR-bar scheme to the
 * on-shell scheme.
 *
 * The function assumes that the physical struct is filled with pole
 * masses and corresponding mixing matrices.  From these quantities,
 * the on-shell model parameters are calculated.
 *
 * @param precision accuracy goal for the conversion
 * @param max_iterations maximum number of iterations
 */
void MSSMNoFV_onshell::convert_to_onshell(
   double precision, unsigned max_iterations)
{
   check_input();
   calculate_DRbar_masses();
   copy_susy_masses_to_pole();

   convert_gauge_couplings();
   convert_BMu();
   convert_vev();
   convert_yukawa_couplings_treelevel();
   convert_Mu_M1_M2(precision, max_iterations);
   convert_yukawa_couplings(); // first guess of resummed yukawas
   convert_ml2();
   convert_yukawa_couplings();
   calculate_MSm(); // calculate MSm with new value of ml2(1,1) and ym
   convert_me2(precision, max_iterations);
   convert_yukawa_couplings();

   // final mass spectrum
   get_problems().clear_problems();
   calculate_DRbar_masses();
   check_problems();
}

/**
 * Calculates the masses.
 *
 * The function assumes that the physical struct is filled with pole
 * masses and corresponding mixing matrices.  From these quantities,
 * the on-shell model parameters are calculated.
 */
void MSSMNoFV_onshell::calculate_masses() {
   check_input();
   convert_gauge_couplings();
   convert_BMu();
   convert_vev();
   convert_yukawa_couplings_treelevel();
   convert_yukawa_couplings(); // tan(beta) resummation in Yukawas

   // final mass spectrum
   get_problems().clear();
   calculate_DRbar_masses();
   check_problems();

   // copy SUSY masses to physical struct
   copy_susy_masses_to_pole();
   get_physical().MAh = get_MAh();
   get_physical().Mhh = get_Mhh();
}

void MSSMNoFV_onshell::check_input() const
{
#define WARN_OR_THROW_IF_ZERO(mass,msg)         \
   if (is_zero(get_##mass())) {                 \
      if (do_force_output())                    \
         WARNING(msg);                          \
      else                                      \
         throw EInvalidInput(msg);              \
   }

   WARN_OR_THROW_IF_ZERO(MW    , "W mass is zero");
   WARN_OR_THROW_IF_ZERO(MZ    , "Z mass is zero");
   WARN_OR_THROW_IF_ZERO(MM    , "Muon mass is zero");
   WARN_OR_THROW_IF_ZERO(Mu    , "mu parameter is zero");
   WARN_OR_THROW_IF_ZERO(MassB , "Bino mass M1 is zero");
   WARN_OR_THROW_IF_ZERO(MassWB, "Wino mass M2 is zero");
   WARN_OR_THROW_IF_ZERO(TB    , "tan(beta) is zero");

#undef WARN_OR_THROW_IF_ZERO
}

void MSSMNoFV_onshell::check_problems() const
{
   if (get_problems().have_problem()) {
      if (!do_force_output())
         throw EPhysicalProblem(get_problems().get_problems());
   }
   if (get_mu2().diagonal().minCoeff() < 0. ||
       get_md2().diagonal().minCoeff() < 0. ||
       get_mq2().diagonal().minCoeff() < 0. ||
       get_me2().diagonal().minCoeff() < 0. ||
       get_ml2().diagonal().minCoeff() < 0.) {
      if (!do_force_output())
         throw EInvalidInput("soft mass squared < 0");
   }
   if (is_zero(get_MCha(0))) {
      if (!do_force_output())
         throw EInvalidInput("lightest chargino mass = 0");
   }
}

void MSSMNoFV_onshell::copy_susy_masses_to_pole()
{
   // if pole masses are zero, interpret the tree-level masses as pole
   // masses

#define COPY_IF_ZERO_0(m)                                               \
   if (MSSMNoFV_onshell::is_zero(get_physical().m))                     \
      get_physical().m = get_##m();
#define COPY_IF_ZERO_1(m,z)                                             \
   if (MSSMNoFV_onshell::is_zero(get_physical().m)) {                   \
      get_physical().m = get_##m();                                     \
      get_physical().z = get_##z();                                     \
   }
#define COPY_IF_ZERO_2(m,u,v)                                           \
   if (MSSMNoFV_onshell::is_zero(get_physical().m)) {                   \
      get_physical().m = get_##m();                                     \
      get_physical().u = get_##u();                                     \
      get_physical().v = get_##v();                                     \
   }

   COPY_IF_ZERO_1(MChi, ZN);
   COPY_IF_ZERO_2(MCha, UM, UP);
   COPY_IF_ZERO_0(MSveL);
   COPY_IF_ZERO_0(MSvmL);
   COPY_IF_ZERO_0(MSvtL);
   COPY_IF_ZERO_1(MSd, ZD);
   COPY_IF_ZERO_1(MSu, ZU);
   COPY_IF_ZERO_1(MSe, ZE);
   COPY_IF_ZERO_1(MSm, ZM);
   COPY_IF_ZERO_1(MStau, ZTau);
   COPY_IF_ZERO_1(MSs, ZS);
   COPY_IF_ZERO_1(MSc, ZC);
   COPY_IF_ZERO_1(MSb, ZB);
   COPY_IF_ZERO_1(MSt, ZT);

#undef COPY_IF_ZERO_0
#undef COPY_IF_ZERO_1
#undef COPY_IF_ZERO_2
}

void MSSMNoFV_onshell::convert_gauge_couplings()
{
   const double MW = get_MW(); // pole mass
   const double MZ = get_MZ(); // pole mass
   const double cW = MW / MZ;  // on-shell weak mixing angle

   set_g1(sqrt(5. / 3.) * EL / cW);
   set_g2(EL / sqrt(1. - sqr(cW)));
}

void MSSMNoFV_onshell::convert_BMu()
{
   const double TB = get_TB(); // DR-bar
   const double MA = get_MA0(); // pole mass
   const double sin2b = 2. * TB / (1. + sqr(TB));

   set_BMu(0.5 * sqr(MA) * sin2b);
}

void MSSMNoFV_onshell::convert_vev()
{
   const double TB = get_TB(); // DR-bar
   const double MW = get_MW(); // pole mass
   const double vev = 2. * MW / get_g2();

   set_vu(vev / sqrt(1. + 1. / sqr(TB)));
   set_vd(get_vu() / TB);
}

/**
 * Returns mb(MZ) in DR-bar scheme.
 *
 * mb(MZ) DR-bar is calculated from mb(mb) MS-bar using the function
 * \a calculate_mb_SM5_DRbar .
 *
 * @return mb(MZ) DR-bar
 */
double MSSMNoFV_onshell::get_MB() const
{
   const double mb_mb = get_MBMB();
   const double alpha_s = calculate_alpha(get_g3());
   const double mb_DRbar = calculate_mb_SM5_DRbar(mb_mb, alpha_s, get_MZ());

   return mb_DRbar;
}

void MSSMNoFV_onshell::convert_yukawa_couplings_treelevel()
{
   Eigen::Matrix<double,3,3> Ye_neu(Eigen::Matrix<double,3,3>::Zero());
   Ye_neu(0, 0) = sqrt(2.) * get_ME() / get_vd();
   Ye_neu(1, 1) = sqrt(2.) * get_MM() / get_vd();
   Ye_neu(2, 2) = sqrt(2.) * get_ML() / get_vd();
   set_Ye(Ye_neu);

   Eigen::Matrix<double,3,3> Yu_neu(Eigen::Matrix<double,3,3>::Zero());
   Yu_neu(0, 0) = sqrt(2.) * get_MU() / get_vu();
   Yu_neu(1, 1) = sqrt(2.) * get_MC() / get_vu();
   Yu_neu(2, 2) = sqrt(2.) * get_MT() / get_vu();
   set_Yu(Yu_neu);

   Eigen::Matrix<double,3,3> Yd_neu(Eigen::Matrix<double,3,3>::Zero());
   Yd_neu(0, 0) = sqrt(2.) * get_MD() / get_vd();
   Yd_neu(1, 1) = sqrt(2.) * get_MS() / get_vd();
   Yd_neu(2, 2) = sqrt(2.) * get_MB() / get_vd();
   set_Yd(Yd_neu);

   // recalculate trilinear couplings with new Yukawas
   set_TYe(Ye_neu * Ae);
   set_TYu(Yu_neu * Au);
   set_TYd(Yd_neu * Ad);
}

/**
 * Calculate Yukawa couplings with resummed tan(beta) corrections
 *
 * @note The resummations depend on the model parameters ml2, me2, Mu,
 * MassB, MassWB, MassG.  Therefore, this routine needs to be called
 * after all these model parameters have been determined.
 */
void MSSMNoFV_onshell::convert_yukawa_couplings()
{
   Eigen::Matrix<double,3,3> Ye_neu(Eigen::Matrix<double,3,3>::Zero());
   Ye_neu(0, 0) = sqrt(2.) * get_ME() / get_vd();
   Ye_neu(1, 1) = sqrt(2.) * get_MM() / get_vd() / (1 + delta_mu_correction(*this));
   Ye_neu(2, 2) = sqrt(2.) * get_ML() / get_vd() / (1 + delta_tau_correction(*this));
   set_Ye(Ye_neu);

   Eigen::Matrix<double,3,3> Yu_neu(Eigen::Matrix<double,3,3>::Zero());
   Yu_neu(0, 0) = sqrt(2.) * get_MU() / get_vu();
   Yu_neu(1, 1) = sqrt(2.) * get_MC() / get_vu();
   Yu_neu(2, 2) = sqrt(2.) * get_MT() / get_vu();
   set_Yu(Yu_neu);

   Eigen::Matrix<double,3,3> Yd_neu(Eigen::Matrix<double,3,3>::Zero());
   Yd_neu(0, 0) = sqrt(2.) * get_MD() / get_vd();
   Yd_neu(1, 1) = sqrt(2.) * get_MS() / get_vd();
   Yd_neu(2, 2) = sqrt(2.) * get_MB() / get_vd() / (1 + delta_bottom_correction(*this));
   set_Yd(Yd_neu);

   // recalculate trilinear couplings with new Yukawas
   set_TYe(Ye_neu * Ae);
   set_TYu(Yu_neu * Au);
   set_TYd(Yd_neu * Ad);
}


template <class Derived>
bool MSSMNoFV_onshell::is_equal(const Eigen::ArrayBase<Derived>& a,
                                const Eigen::ArrayBase<Derived>& b,
                                double precision_goal)
{
   return (a - b).cwiseAbs().maxCoeff() < precision_goal;
}

bool MSSMNoFV_onshell::is_equal(double a, double b, double precision_goal)
{
   return flexiblesusy::is_equal(a, b, precision_goal);
}

template <class Derived>
bool MSSMNoFV_onshell::is_zero(const Eigen::ArrayBase<Derived>& a,
                               double eps)
{
   return a.cwiseAbs().maxCoeff() < eps;
}

bool MSSMNoFV_onshell::is_zero(double a, double eps)
{
   return flexiblesusy::is_zero(a, eps);
}

/**
 * Returns index of most bino-like neutralino.  The function extracts
 * this information from the given neutralino mixing matrix.
 *
 * @param ZN neutralino mixing matrix
 */
template <class Derived>
unsigned MSSMNoFV_onshell::find_bino_like_neutralino(
   const Eigen::MatrixBase<Derived>& ZN)
{
   unsigned max_bino;
   ZN.col(0).cwiseAbs().maxCoeff(&max_bino);

   return max_bino;
}

/**
 * Returns index of most right-handed smuon.  The function extracts
 * this information from the given smuon mixing matrix.
 *
 * @param ZM smuon mixing matrix
 */
template <class Derived>
unsigned MSSMNoFV_onshell::find_right_like_smuon(
   const Eigen::MatrixBase<Derived>& ZM)
{
   return (ZM(0,0) > ZM(0,1)) ? 1 : 0;
}

/**
 * Determines the Mu parameter and the soft-breaking Bino and Wino
 * mass parameters from the two chargino pole masses and the most
 * bino-like neutralino pole mass.  The function uses a fixed-point
 * iteration.
 *
 * @param precision_goal precision goal of iteration
 * @param max_iterations maximum number of iterations
 */
void MSSMNoFV_onshell::convert_Mu_M1_M2(
   double precision_goal,
   unsigned max_iterations)
{
   // find neutralino, which is most bino like
   const unsigned max_bino = find_bino_like_neutralino(get_physical().ZN);

   const auto MCha_goal(get_physical().MCha);
   auto MChi_goal(get_MChi());
   MChi_goal(max_bino) = get_physical().MChi(max_bino);

   if (verbose_output) {
      std::cout << "Converting Mu, M1, M2 to on-shell scheme ...\n"
                   "   Goal: MCha = " << MCha_goal.transpose()
                << ", MChi(" << max_bino << ") = " << MChi_goal(max_bino)
                << '\n';
   }

   bool accuracy_goal_reached =
      MSSMNoFV_onshell::is_equal(MCha_goal, get_MCha(), precision_goal) &&
      MSSMNoFV_onshell::is_equal(MChi_goal(max_bino), get_MChi(max_bino), precision_goal);
   unsigned it = 0;

   while (!accuracy_goal_reached && it < max_iterations) {

      const auto U(get_UM()); // neg. chargino mixing matrix
      const auto V(get_UP()); // pos. chargino mixing matrix
      const auto N(get_ZN()); // neutralino mixing matrix
      const auto X(U.transpose() * MCha_goal.matrix().asDiagonal() * V);
      const auto Y(N.transpose() * MChi_goal.matrix().asDiagonal() * N);

      set_MassB(std::real(Y(0,0)));
      set_MassWB(std::real(X(0,0)));
      set_Mu(std::real(X(1,1)));

      calculate_DRbar_masses();

      MChi_goal = get_MChi();
      MChi_goal(max_bino) = get_physical().MChi(max_bino);

      if (verbose_output) {
         std::cout << "   Iteration " << it << ": Mu = " << get_Mu()
                   << ", M1 = " << get_MassB()
                   << ", M2 = " << get_MassWB()
                   << ", MCha = " << get_MCha().transpose()
                   << ", MChi(" << max_bino << ") = " << get_MChi(max_bino)
                   << '\n';
      }

      accuracy_goal_reached =
         MSSMNoFV_onshell::is_equal(MCha_goal, get_MCha(), precision_goal) &&
         MSSMNoFV_onshell::is_equal(MChi_goal(max_bino), get_MChi(max_bino), precision_goal);

      it++;
   }

   if (it == max_iterations) {
      const double precision =
         std::max((MCha_goal - get_MCha()).cwiseAbs().maxCoeff(),
                  std::abs(MChi_goal(max_bino) - get_MChi(max_bino)));
      get_problems().flag_no_convergence_Mu_MassB_MassWB(precision, it);
   } else {
      get_problems().unflag_no_convergence_Mu_MassB_MassWB();
   }

   if (verbose_output) {
      const double precision =
         std::max((MCha_goal - get_MCha()).cwiseAbs().maxCoeff(),
                  std::abs(MChi_goal(max_bino) - get_MChi(max_bino)));
      if (it == max_iterations) {
         std::cout <<
            "   DR-bar to on-shell conversion for Mu, M1 and M2 did"
            " not converge (reached absolute accuracy: " << precision <<
            " GeV, accuracy goal: " << precision_goal <<
            ", max. iterations: " << max_iterations << ")\n";
      }
      std::cout << "   Achieved absolute accuracy: " << precision << " GeV\n";
   }
}

/**
 * Determines soft-breaking left-handed smuon mass parameter from the
 * muon sneutrino pole mass.
 */
void MSSMNoFV_onshell::convert_ml2()
{
   const double MSvmL_pole = get_physical().MSvmL;
   const double vd2 = sqr(get_vd());
   const double vu2 = sqr(get_vu());
   const double g12 = sqr(get_g1());
   const double g22 = sqr(get_g2());

   if (verbose_output) {
      std::cout << "Converting msl(2,2) to on-shell scheme ...\n"
         "   Using MSvm_pole = " << MSvmL_pole << '\n';
   }

   // calculate ml2(1,1) from muon sneutrino pole mass
   const double ml211
      = sqr(MSvmL_pole) + 0.125*(0.6*g12*(vu2 - vd2) + g22*(vu2 - vd2));

   set_ml2(1,1,ml211);
   calculate_MSvmL();

   if (verbose_output) {
      std::cout << "   New msl(2,2) = " << ml211
                << ", new MSvmL = " << get_MSvmL() << '\n';
   }
}

/**
 * Determines soft-breaking right-handed smuon mass parameter from one
 * smuon pole mass.
 *
 * @param precision_goal precision goal of iteration
 * @param max_iterations maximum number of iterations
 */
void MSSMNoFV_onshell::convert_me2(
   double precision_goal,
   unsigned max_iterations)
{
   double precision = convert_me2_fpi(precision_goal, max_iterations);

   if (precision > precision_goal)
      precision = convert_me2_root(precision_goal, max_iterations);

   if (precision > precision_goal)
      get_problems().flag_no_convergence_me2(precision, max_iterations);
   else
      get_problems().unflag_no_convergence_me2();
}

/**
 * Determines soft-breaking right-handed smuon mass parameter from one
 * smuon pole mass.  The function uses a one-dimensional root-finding
 * algorithm.
 *
 * @param precision_goal precision goal for the root finding algorithm
 * @param max_iterations maximum number of iterations
 * @return achieved precision
 */
double MSSMNoFV_onshell::convert_me2_root_modify(
   double precision_goal,
   unsigned max_iterations)
{
   class Difference_MSm {
   public:
      Difference_MSm(const MSSMNoFV_onshell& model_)
         : model(model_) {}

      double operator()(double me211) {
         model.set_me2(1,1,me211);
         model.calculate_MSm();

         Eigen::Array<double,2,1> MSm_pole(model.get_physical().MSm);
         std::sort(MSm_pole.data(), MSm_pole.data() + MSm_pole.size());
         const int right_index = find_right_like_smuon(model.get_ZM());

         return model.get_MSm(right_index) - MSm_pole(right_index);
      }
   private:
      MSSMNoFV_onshell model;
   };

   boost::uintmax_t it = max_iterations;

   // stopping criterion, given two brackets a, b
   auto Stop_crit = [precision_goal](double a, double b) -> bool {
      return flexiblesusy::is_equal(a,b,precision_goal);
   };

   if (verbose_output) {
      std::cout << "Converting mse(2,2) to on-shell scheme with "
         "root finder ...\n";
   }

   // initial guess for the brackets of me2(1,1)
   const double initial_bracket = sqr(1e3 * get_MSm().cwiseAbs().maxCoeff());

   // find the root
   try {
      const std::pair<double,double> root = boost::math::tools::toms748_solve(
         Difference_MSm(*this), 0., initial_bracket, Stop_crit, it);
      set_me2(1, 1, 0.5*(root.first + root.second));
   } catch (const std::exception& e) {
      if (verbose_output) {
         std::cout <<
            "   DR-bar to on-shell conversion for mse failed with"
            " root finder: " << e.what() << '\n';
      }
   }

   calculate_MSm();

   const double precision = std::abs(Difference_MSm(*this)(get_me2(1,1)));

   if (verbose_output) {
      if (it >= max_iterations) {
         std::cout <<
            "   DR-bar to on-shell conversion for mse did not converge with"
            " root finder (reached absolute accuracy: " << precision <<
            " GeV, accuracy goal: " << precision_goal <<
            ", max. iterations: " << max_iterations << ")\n";
      }
      std::cout << "   Achieved absolute accuracy: "
                << precision << " GeV\n";
   }

   return precision;
}

/**
 * Determines soft-breaking right-handed smuon mass parameter from one
 * smuon pole mass.  The function uses a one-dimensional root-finding
 * algorithm.
 *
 * If a NaN appears during the FPI, the function resets the
 * soft-breaking right-handed smuon mass parameter to the initial
 * value and returns the maximum double.
 *
 * @param precision_goal precision goal for the root finding algorithm
 * @param max_iterations maximum number of iterations
 * @return achieved precision
 */
double MSSMNoFV_onshell::convert_me2_root(
   double precision_goal,
   unsigned max_iterations)
{
   const double me2_save = get_me2(1,1);
   const double precision = convert_me2_root_modify(precision_goal, max_iterations);

   if (!std::isfinite(precision) ||
       !std::isfinite(get_me2(1,1)) ||
       !get_MSm().allFinite() ||
       !get_ZM().allFinite()) {
      // reset
      set_me2(1,1,me2_save);
      calculate_MSm();

      return std::numeric_limits<double>::max();
   }

   return precision;
}

/**
 * Determines soft-breaking right-handed smuon mass parameter from one
 * smuon pole mass.  The function uses a fixed-point iteration.
 *
 * @param precision_goal precision goal of iteration
 * @param max_iterations maximum number of iterations
 * @return achieved precision
 */
double MSSMNoFV_onshell::convert_me2_fpi_modify(
   double precision_goal,
   unsigned max_iterations)
{
   // sorted pole masses
   Eigen::Array<double,2,1> MSm_pole_sorted(get_physical().MSm);
   std::sort(MSm_pole_sorted.data(),
             MSm_pole_sorted.data() + MSm_pole_sorted.size());

   Eigen::Array<double,2,1> MSm_goal(MSm_pole_sorted);

   int right_index = find_right_like_smuon(get_ZM());

   if (verbose_output) {
      std::cout << "Converting mse(2,2) to on-shell scheme with FPI ...\n"
                   "   Goal: MSm(" << right_index << ") = "
                << MSm_goal(right_index) << '\n';
   }

   bool accuracy_goal_reached =
      MSSMNoFV_onshell::is_equal(get_MSm(right_index), MSm_goal(right_index),
                                 precision_goal);
   unsigned it = 0;

   while (!accuracy_goal_reached && it < max_iterations) {
      const Eigen::Matrix<double,2,2> ZM(get_ZM()); // smuon mixing matrix
      const Eigen::Matrix<double,2,2> M(
         ZM.adjoint() * MSm_goal.square().matrix().asDiagonal() * ZM);

      const double vd2 = sqr(get_vd());
      const double vu2 = sqr(get_vu());
      const double g12 = sqr(get_g1());
      const double ymu2 = std::norm(Ye(1,1));

      const double me211 = M(1,1)
         - (0.5*ymu2*vd2 - 0.15*g12*vd2 + 0.15*g12*vu2);

      set_me2(1,1,me211);
      calculate_MSm();

      if (!std::isfinite(me211) || !get_MSm().allFinite() || !get_ZM().allFinite()) {
         if (verbose_output) {
            std::cout << "   NaN appearing in DR-bar to on-shell conversion"
               " for mse with FPI:\n"
               "      mse2(2,2) = " << me211 <<
               ", MSm = " << get_MSm().transpose() <<
               ", ZM = " << get_ZM().row(0) << ' ' << get_ZM().row(1) << '\n';
         }
         return std::numeric_limits<double>::max();
      }

      right_index = find_right_like_smuon(get_ZM());

      MSm_goal = get_MSm();
      MSm_goal(right_index) = MSm_pole_sorted(right_index);

      if (verbose_output) {
         std::cout << "   Iteration " << it << ": mse(2,2) = "
                   << signed_abs_sqrt(me211) << ", MSm(" << right_index
                   << ") = " << get_MSm(right_index) << '\n';
      }

      accuracy_goal_reached =
         MSSMNoFV_onshell::is_equal(get_MSm(right_index), MSm_goal(right_index),
                                    precision_goal);

      it++;
   }

   const double precision =
      std::abs(get_MSm(right_index) - MSm_goal(right_index));

   if (verbose_output) {
      if (it == max_iterations) {
         std::cout <<
            "   DR-bar to on-shell conversion for mse did not converge with"
            " FPI (reached absolute accuracy: " << precision <<
            " GeV, accuracy goal: " << precision_goal <<
            ", max. iterations: " << max_iterations << ")\n";
      }
      std::cout << "   Achieved absolute accuracy: "
                << precision << " GeV\n";
   }

   return precision;
}

/**
 * Determines soft-breaking right-handed smuon mass parameter from one
 * smuon pole mass.  The function uses a fixed-point iteration.
 *
 * If a NaN appears during the FPI, the function resets the
 * soft-breaking right-handed smuon mass parameter to the initial
 * value and returns the maximum double.
 *
 * @param precision_goal precision goal of iteration
 * @param max_iterations maximum number of iterations
 * @return achieved precision
 */
double MSSMNoFV_onshell::convert_me2_fpi(
   double precision_goal,
   unsigned max_iterations)
{
   const double me2_save = get_me2(1,1);
   const double precision = convert_me2_fpi_modify(precision_goal, max_iterations);

   if (!std::isfinite(precision) ||
       !std::isfinite(get_me2(1,1)) ||
       !get_MSm().allFinite() ||
       !get_ZM().allFinite()) {
      // reset
      set_me2(1,1,me2_save);
      calculate_MSm();

      return std::numeric_limits<double>::max();
   }

   return precision;
}

std::ostream& operator<<(std::ostream& os, const MSSMNoFV_onshell& model)
{
   os <<
      "======================================\n"
      " (g-2) parameters \n"
      "======================================\n"
      << "1/alpha(MZ) = " << 1./calculate_alpha(model.get_EL()) << '\n'
      << "1/alpha(0)  = " << 1./calculate_alpha(model.get_EL0()) << '\n'
      << "alpha_s(MZ) = " <<    calculate_alpha(model.get_g3()) << '\n'
      <<
      "--------------------------------------\n"
      " on-shell masses and parameters \n"
      "--------------------------------------\n"
      "MM          = " << model.get_MM() << '\n' <<
      "MT          = " << model.get_MT() << '\n' <<
      "mb(mb)      = " << model.get_MBMB() << '\n' <<
      "mb(MZ)      = " << model.get_MB() << '\n' <<
      "MTau        = " << model.get_ML() << '\n' <<
      "MW          = " << model.get_MW() << '\n' <<
      "MZ          = " << model.get_MZ() << '\n' <<
      "MSm         = " << model.get_MSm().transpose() << '\n' <<
      "USm         = " << model.get_USm().row(0) << ' '
                       << model.get_USm().row(1) << '\n' <<
      "MSvm        = " << model.get_MSvmL() << '\n' <<
      "MSb         = " << model.get_MSb().transpose() << '\n' <<
      "USb         = " << model.get_USb().row(0) << ' '
                       << model.get_USb().row(1) << '\n' <<
      "MSt         = " << model.get_MSt().transpose() << '\n' <<
      "USt         = " << model.get_USt().row(0) << ' '
                       << model.get_USt().row(1) << '\n' <<
      "MStau       = " << model.get_MStau().transpose() << '\n' <<
      "UStau       = " << model.get_UStau().row(0) << ' '
                       << model.get_UStau().row(1) << '\n' <<
      "MCha        = " << model.get_MCha().transpose() << '\n' <<
      "UM          = " << model.get_UM().row(0) << ' '
                       << model.get_UM().row(1) << '\n' <<
      "UP          = " << model.get_UP().row(0) << ' '
                       << model.get_UP().row(1) << '\n' <<
      "MChi        = " << model.get_MChi().transpose() << '\n' <<
      "ZN          = " << model.get_ZN().row(0) << '\n' <<
      "              " << model.get_ZN().row(1) << '\n' <<
      "              " << model.get_ZN().row(2) << '\n' <<
      "              " << model.get_ZN().row(3) << '\n' <<
      "MA0         = " << model.get_MA0() << '\n' <<
      "MH          = " << model.get_Mhh().transpose() << '\n' <<
      "tan(beta)   = " << model.get_TB() << '\n' <<
      "yu          = " << model.get_Yu().diagonal().transpose() << '\n' <<
      "yd resummed = " << model.get_Yd().diagonal().transpose() << '\n' <<
      "ye resummed = " << model.get_Ye().diagonal().transpose() << '\n' <<
      "yd non res. = " << (sqrt(2.) * model.get_MD() / model.get_vd()) <<
                   " " << (sqrt(2.) * model.get_MS() / model.get_vd()) <<
                   " " << (sqrt(2.) * model.get_MB() / model.get_vd()) << '\n' <<
      "ye non res. = " << (sqrt(2.) * model.get_ME() / model.get_vd()) <<
                   " " << (sqrt(2.) * model.get_MM() / model.get_vd()) <<
                   " " << (sqrt(2.) * model.get_ML() / model.get_vd()) << '\n' <<
      "Mu          = " << model.get_Mu() << '\n' <<
      "M1          = " << model.get_MassB() << '\n' <<
      "M2          = " << model.get_MassWB() << '\n' <<
      "M3          = " << model.get_MassG() << '\n' <<
      "msl         = " << model.get_ml2().diagonal().transpose().unaryExpr(std::ptr_fun(signed_abs_sqrt)) << '\n' <<
      "mse         = " << model.get_me2().diagonal().transpose().unaryExpr(std::ptr_fun(signed_abs_sqrt)) << '\n' <<
      "msq         = " << model.get_mq2().diagonal().transpose().unaryExpr(std::ptr_fun(signed_abs_sqrt)) << '\n' <<
      "msu         = " << model.get_mu2().diagonal().transpose().unaryExpr(std::ptr_fun(signed_abs_sqrt)) << '\n' <<
      "msd         = " << model.get_md2().diagonal().transpose().unaryExpr(std::ptr_fun(signed_abs_sqrt)) << '\n' <<
      "Au          = " << model.get_Au().diagonal().transpose() << '\n' <<
      "Ad          = " << model.get_Ad().diagonal().transpose() << '\n' <<
      "Ae          = " << model.get_Ae().diagonal().transpose() << '\n' <<
      "ren. scale  = " << model.get_scale() << '\n'
      ;

   return os;
}

} // gm2calc
