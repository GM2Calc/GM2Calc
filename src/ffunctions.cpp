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

#include "ffunctions.hpp"
#include "dilog.hpp"
#include "numerics2.hpp"

#include <iostream>
#include <cmath>
#include <complex>

#define ERROR(message) std::cerr << "Error: " << message << '\n';

namespace gm2calc {

namespace {
   const double Pi = 3.14159265358979323846;

   /// returns number cubed
   template <typename T> T pow3(T x) { return x*x*x; }

   /// returns number to the power 4
   template <typename T> T pow4(T x) { return x*x*x*x; }
}

using namespace flexiblesusy;

double abs_sqrt(double x) {
   return std::sqrt(std::abs(x));
}

int sign(double x) { return x < 0 ? -1 : 1; }

double signed_sqr(double x) { return sign(x) * x * x; }

double signed_abs_sqrt(double x) {
   return sign(x) * std::sqrt(std::abs(x));
}

double F1C(double x) {
   if (is_zero(x))
      return 4.;

   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.05)) {
      return 1.0 + d*(-0.6 + d*(0.4 + d*(-2.0/7.0
         + d*(3.0/14.0 + d*(-1.0/6.0
         + 2.0/15.0*d)))));
   }

   return 2. / pow4(d) * (2. + 3. * x - 6. * sqr(x)
                          + pow3(x) + 6. * x * std::log(x));
}

double F2C(double x) {
   if (is_equal(x, 1., 0.03)) {
      return 1. - 0.75*(-1 + x) + 0.6*sqr(-1 + x) - 0.5*pow3(-1 + x)
         + 3/7.*pow4(-1 + x) - 0.375*std::pow(-1 + x,5)
         + 1/3.*std::pow(-1 + x,6);
   }

   if (is_zero(x))
      return 0.;

   return 3. / (2. * pow3(1. - x)) * (- 3. + 4. * x - sqr(x) - 2. * log(x));
}

double F3C(double x) {
   if (is_equal(x, 1., 0.07)) {
      return 1. + 1059/1175.*(-1 + x)
         - 4313/3525.*sqr(-1 + x)
         + 70701/57575.*pow3(-1 + x)
         - 265541/230300.*pow4(-1 + x)
         + 48919/46060.*std::pow(-1 + x,5)
         - 80755/82908.*std::pow(-1 + x,6);
   }

   return ( 4. / (141. * pow4(1. - x)) * ((1. - x) * (151. * sqr(x) - 335. * x + 592.)
             + 6. * (21. * pow3(x) - 108. * sqr(x) - 93. * x + 50.) * log(x)
             - 54. * x * (sqr(x) - 2. * x - 2.) * sqr(log(x))
             - 108. * x * (sqr(x) - 2. * x + 12.) * dilog(1.- x)) );
}

double F4C(double x) {
   if (is_equal(x, 1., 0.07)) {
      return 1. - 45/122.*(-1 + x)
         + 941/6100.*sqr(-1 + x)
         - 17/305.*pow3(-1 + x)
         + 282/74725.*pow4(-1 + x)
         + 177/6832.*std::pow(-1 + x,5)
         - 47021/1076040.*std::pow(-1 + x,6);
   }

   if (is_zero(x))
      return 0.;

   return ( - 9. / (122. * pow3(1. - x)) * (8. * (sqr(x) - 3. * x + 2.)
             +  (11. * sqr(x) - 40. * x + 5.) * log(x)
             - 2. * (sqr(x) - 2. * x - 2.) * sqr(log(x))
             - 4. * (sqr(x) - 2. * x + 9.) * dilog(1.- x)) );
}

double F1N(double x) {
   if (is_equal(x, 1., 0.06)) {
      return 1. - 0.4*(-1 + x) + 0.2*sqr(-1 + x) - 4/35.*pow3(-1 + x)
         + 1/14.*pow4(-1 + x) - 1/21.*std::pow(-1 + x,5)
         + 1/30.*std::pow(-1 + x,6);
   }

   if (is_zero(x))
      return 2.;

   return 2. / pow4(1. - x) * (1. - 6. * x + 3. * sqr(x)
                               + 2. * pow3(x) - 6. * sqr(x) * log(x));
}

double F2N(double x) {
   if (is_equal(x, 1., 0.04)) {
      return 1. - 0.5*(-1 + x) + 0.3*sqr(-1 + x) - 0.2*pow3(-1 + x)
         + 1/7.*pow4(-1 + x) - 3/28.*std::pow(-1 + x,5)
         + 1/12.*std::pow(-1 + x,6);
   }

   if (is_zero(x))
      return 3.;

   return 3. / pow3(1. - x) * (1. - sqr(x) + 2. * x * log(x));
}

double F3N(double x) {
   if (is_equal(x, 1., 0.07)) {
      return 1. + 76/875.*(-1 + x)
         - 431/2625.*sqr(-1 + x)
         + 5858/42875.*pow3(-1 + x)
         - 3561/34300.*pow4(-1 + x)
         + 23/294.*std::pow(-1 + x,5)
         - 4381/73500.*std::pow(-1 + x,6);
   }

   if (is_zero(x))
      return 8./105.;

   return 4. / (105. * pow4(1. - x)) * ((1. - x) * (- 97. * sqr(x) - 529. * x + 2.)
            + 6. * sqr(x) * (13. * x + 81.) * log(x)
            + 108. * x * (7. * x + 4.) * dilog(1. - x));
}

double F4N(double x) {
   if (is_equal(x, 1., 0.07)) {
      return 1. - 0.13875*sqr(-1 + x)
         + 0.1475*pow3(-1 + x)
         - 129/980.*pow4(-1 + x)
         + 177/1568.*std::pow(-1 + x,5)
         - 775/8064.*std::pow(-1 + x,6);
   }

   if (is_zero(x))
      return -(3./4.)*(-9. + sqr(Pi));

   return - 2.25 / pow3(1. - x) * ((x + 3.) * (x * log(x) + x - 1.)
                                   + (6. * x + 2.) * dilog(1. - x));
}

/// Fb(1,1)
double Fb11(double x, double y) {
   return (293. - 216.*y + 63.*sqr(y) + sqr(x)*(63. - 50.*y + 15.*sqr(y))
           - 2.*x*(108. - 84.*y + 25.*sqr(y))) / 840.;
}

/// Fb(x,1)
double Fb1(double x, double y) {
   return (2. + 3.*x - 6.*sqr(x) + pow3(x) + 6.*x*log(x)) / (6.*pow4(-1. + x))
      + (-1. + y) * (3. + 10.*x - 18.*sqr(x) + 6.*pow3(x) - pow4(x) + 12.*x*log(x))
      / (12.*std::pow(-1. + x,5))
      + sqr(-1. + y) * (12. + 65.*x - 120.*sqr(x) + 60.*pow3(x) - 20.*pow4(x)
                        + 3.*std::pow(x,5) + 60.*x*log(x)) / (60.*std::pow(-1. + x,6));
}

/// Fb(x,x)
double Fbx(double x, double y) {
   return (-5. + 4.*x + sqr(x) - 2.*log(x) - 4.*x*log(x)) / (2.*pow4(-1. + x))
      - (-x + y) * (-1. - 9.*x + 9.*sqr(x) + pow3(x) - 6.*x*log(x) - 6.*sqr(x)*log(x))
      / (2.*std::pow(-1. + x,5)*x)
      - sqr(-x + y) * (-1. + 12.*x + 36.*sqr(x) - 44.*pow3(x) - 3.*pow4(x) + 36.*sqr(x)*log(x)
                       + 24.*pow3(x)*log(x)) / (6.*std::pow(-1. + x,6)*sqr(x));
}

double Fb(double x, double y) {
   if ((is_zero(x) && is_zero(y)) || is_zero(x) || is_zero(y))
      return 0.;

   if (is_equal(x, 1., 0.01) && is_equal(y, 1., 0.01))
      return Fb11(x,y);

   if (is_equal(x, 1., 0.01))
      return Fb1(y,x);

   if (is_equal(y, 1., 0.01))
      return Fb1(x,y);

   if (is_equal_rel(x, y, 0.01))
      return Fbx(x,y);

   return - (G4(x) - G4(y)) / (x - y);
}

/// Fa(1,1)
double Fa11(double x, double y) {
   return (1311. - 1158.*y + 365.*sqr(y) + 5.*sqr(x)*(73. - 66.*y + 21.*sqr(y))
           - 2.*x*(579. - 520.*y + 165.*sqr(y))) / 840.;
}

/// Fa(x,1)
double Fa1(double x, double y) {
   return (-11. + 18.*x - 9.*sqr(x) + 2.*pow3(x) - 6.*log(x)) / (6.*pow4(-1. + x))
      + (-1. + y) * (-25. + 48.*x - 36.*sqr(x) + 16.*pow3(x) - 3.*pow4(x) - 12.*log(x))
      / (12.*std::pow(-1. + x,5))
      + sqr(-1. + y) * (-137. + 300.*x - 300.*sqr(x) + 200.*pow3(x) - 75.*pow4(x)
                        + 12.*std::pow(x,5) - 60.*log(x)) / (60.*std::pow(-1. + x,6));
}

/// Fa(x,x)
double Fax(double x, double y) {
   return (2. + 3.*x - 6.*sqr(x) + pow3(x) + 6.*x*log(x)) / (2.*pow4(-1. + x)*x)
      - (-x + y) * (-1. + 8.*x - 8.*pow3(x) + pow4(x) + 12.*sqr(x)*log(x))
      / (2.*std::pow(-1. + x,5)*sqr(x))
      - sqr(-x + y) * (-2. + 15.*x - 60.*sqr(x) + 20.*pow3(x) + 30.*pow4(x) - 3.*std::pow(x,5)
                       - 60.*pow3(x)*log(x)) / (6.*std::pow(-1. + x,6)*pow3(x));
}

double Fa(double x, double y) {
   if ((is_zero(x) && is_zero(y)) || is_zero(x) || is_zero(y))
      return 0.;

   if (is_equal(x, 1., 0.01) && is_equal(y, 1., 0.01))
      return Fa11(x,y);

   if (is_equal(x, 1., 0.01))
      return Fa1(y,x);

   if (is_equal(y, 1., 0.01))
      return Fa1(x,y);

   if (is_equal_rel(x, y, 0.01))
      return Fax(x,y);

   return - (G3(x) - G3(y)) / (x - y);
}

double G3(double x) {
   return 1. / (2. * pow3(x - 1.)) * ((x - 1.) * (x - 3.) + 2. * log(x));
}

double G4(double x) {
   return 1. / (2. * pow3(x - 1.)) * ((x - 1.) * (x + 1.) - 2. * x * log(x));
}

/// Iabc(a,a,a)
double Iaaa(double a, double b, double c) {
   return (151.*pow4(a) + 13.*sqr(b)*sqr(c) - 128.*pow3(a)*(b + c) - 40.*a*b*c*(b + c)
           + sqr(a)*(37.*sqr(b) + 128.*b*c + 37.*sqr(c))) / (60.*std::pow(a,6));
}

/// Iabc(a,a,c)
double Iaac(double a, double b, double c) {
   return ((sqr(a) - sqr(c))
           * (17.*std::pow(a,6) - 16.*std::pow(a,5)*b - 40.*pow3(a)*b*sqr(c)
              + 8.*a*b*pow4(c) - sqr(b)*pow4(c) + pow4(a)*(5.*sqr(b) + 8.*sqr(c))
              + sqr(a)*(20.*sqr(b)*sqr(c) - pow4(c)))
           - 6.*sqr(a)*sqr(c) * log(sqr(a)/sqr(c))
           * (6.*pow4(a) - 8.*pow3(a)*b + 3.*sqr(a)*(sqr(b) - sqr(c)) + sqr(c)*(sqr(b) + sqr(c))))
      / (6.*sqr(a)*pow4(sqr(a) - sqr(c)));
}

/// Iabc(a,a,0)
double Iaa0(double a, double b) {
   return (17.*sqr(a) - 16.*a*b + 5.*sqr(b)) / (6.*pow4(a));
}

/// Iabc(0,b,c)
double I0bc(double b, double c) {
   return log(sqr(b/c))/(sqr(b) - sqr(c));
}

double Iabc(double a, double b, double c) {
   if ((is_zero(a) && is_zero(b) && is_zero(c)) ||
       (is_zero(a) && is_zero(b)) ||
       (is_zero(a) && is_zero(c)) ||
       (is_zero(b) && is_zero(c)))
      return 0.;

   if (is_equal_rel(std::abs(a), std::abs(b), 0.01) && is_equal_rel(std::abs(a), std::abs(c), 0.01))
      return Iaaa(std::abs(a),std::abs(b),std::abs(c));

   if (is_equal_rel(std::abs(a), std::abs(b), 0.01)) {
      if (is_zero(c))
         return Iaa0(std::abs(a),std::abs(b));
      return Iaac(std::abs(a),std::abs(b),c);
   }

   if (is_equal_rel(std::abs(b), std::abs(c), 0.01)) {
      if (is_zero(a))
         return Iaa0(std::abs(b),std::abs(c));
      return Iaac(std::abs(b),std::abs(c),a);
   }

   if (is_equal_rel(std::abs(a), std::abs(c), 0.01)) {
      if (is_zero(b))
         return Iaa0(std::abs(a),std::abs(c));
      return Iaac(std::abs(a),std::abs(c),b);
   }

   if (is_zero(a))
      return I0bc(b,c);

   if (is_zero(b))
      return I0bc(c,a);

   if (is_zero(c))
      return I0bc(a,b);

   return ( (sqr(a * b) * log(sqr(a / b))
           + sqr(b * c) * log(sqr(b / c))
           + sqr(c * a) * log(sqr(c / a)))
           / ((sqr(a) - sqr(b)) * (sqr(b) - sqr(c)) * (sqr(a) - sqr(c))) );
}

/**
 * Calculates \f$f_{PS}(z)\f$, Eq (70) arXiv:hep-ph/0609168
 */
double f_PS(double z) {
   double result = 0.;

   if (z < 0.25) {
      const double y = sqrt(1. - 4. * z);
      result = 2. * z / y * (dilog(1. - 0.5 * (1. - y) / z) - dilog(1. - 0.5 * (1. + y) / z));
   } else {
      const std::complex<double> y(sqrt(std::complex<double>(1. - 4. * z, 0.)));
      const std::complex<double> zc(z, 0.);
      result = std::real(2. * zc / y * (dilog(1. - 0.5 * (1. - y) / zc) - dilog(1. - 0.5 * (1. + y) / zc)));
   }

   return result;
}

/**
 * Calculates \f$f_S(z)\f$, Eq (71) arXiv:hep-ph/0609168
 */
double f_S(double z) {
   if(z < 0.)
      ERROR("f_S: z must not be negativ!");

   return (2. * z - 1.) * f_PS(z) - 2. * z * (2. + log(z));
}

/**
 * Calculates \f$f_{\tilde{f}}(z)\f$, Eq (72) arXiv:hep-ph/0609168
 */
double f_sferm(double z) {
   if(z < 0.)
      ERROR("f_sferm: z must not be negativ!");

   return 0.5 * z * (2. + log(z) - f_PS(z));
}

} // namespace gm2calc
