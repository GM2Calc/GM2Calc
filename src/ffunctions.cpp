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

#include <cmath>
#include <complex>
#include <iostream>

#define ERROR(message) std::cerr << "Error: " << message << '\n';

namespace gm2calc {

namespace {
   /// returns number cubed
   template <typename T> T pow3(T x) { return x*x*x; }

   /// returns number to the power 4
   template <typename T> T pow4(T x) { return x*x*x*x; }
} // anonymous namespace

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
   if (is_zero(x)) {
      return 4.0;
   }

   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.05)) {
      return 1.0 + d*(-0.6 + d*(0.4 + d*(-2.0/7.0
         + d*(3.0/14.0 + d*(-1.0/6.0
         + 2.0/15.0*d)))));
   }

   return 2.0/pow4(d)*(2.0 + x*(3.0 + 6.0*std::log(x) + x*(-6.0 + x)));
}

double F2C(double x) {
   if (is_zero(x)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.03)) {
      const double d = x - 1.0;

      return 1.0 + d*(-0.75 + d*(0.6 + d*(-0.5 + d*(3.0/7.0
         + d*(-0.375 + 1.0/3.0*d)))));
   }

   return 3.0/(2.0*pow3(1.0 - x))*(-3.0 - 2.0*std::log(x) + x*(4.0 - x));
}

double F3C(double x) {
   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.07)) {
      return 1.0
         + d*(1059.0/1175.0
         + d*(-4313.0/3525.0
         + d*(70701.0/57575.0
         + d*(-265541.0/230300.0
         + d*(+48919.0/46060.0
         - 80755.0/82908.0*d)))));
   }

   const double lx = std::log(x);

   return 4.0/(141.0*pow4(d)) * (
      + (1.0 - x) * (151.0 * sqr(x) - 335.0 * x + 592.0)
      + 6.0 * (21.0 * pow3(x) - 108.0 * sqr(x) - 93.0 * x + 50.0) * lx
      - 54.0 * x * (sqr(x) - 2.0 * x - 2.0) * sqr(lx)
      - 108.0 * x * (sqr(x) - 2.0 * x + 12.0) * dilog(1.0 - x)
      );
}

double F4C(double x) {
   if (is_zero(x)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.07)) {
      const double d = x - 1.0;

      return 1.0
         + d*(-45.0/122.0
         + d*(941.0/6100.0
         + d*(-17.0/305.0
         + d*(+282.0/74725.0
         + d*(+177.0/6832.0 - 47021.0/1076040.0*d)))));
   }

   const double lx = std::log(x);

   return -9.0/(122.0 * pow3(1.0 - x)) * (
      + 8.0 * (sqr(x) - 3.0 * x + 2.0)
      + (11.0 * sqr(x) - 40.0 * x + 5.0) * lx
      - 2.0 * (sqr(x) - 2.0 * x - 2.0) * sqr(lx)
      - 4.0 * (sqr(x) - 2.0 * x + 9.0) * dilog(1.0- x)
      );
}

double F1N(double x) {
   if (is_zero(x)) {
      return 2.0;
   }

   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.06)) {
      return 1.0 + d*(-0.4 + d*(0.2 + d*(-4.0/35.0
         + d*(1.0/14.0 + d*(-1.0/21.0 + 1.0/30.0*d)))));
   }

   return 2.0/pow4(d)*(1.0 + x*(-6.0 + x*(+3.0 - 6.0 * std::log(x) + 2.0 * x)));
}

double F2N(double x) {
   if (is_zero(x)) {
      return 3.0;
   }

   if (is_equal(x, 1.0, 0.04)) {
      const double d = x - 1.0;

      return 1. + d*(-0.5 + d*(0.3 + d*(-0.2
         + d*(1.0/7.0 + d*(-3.0/28.0 + 1.0/12.0*d)))));
   }

   return 3.0/pow3(1.0 - x) * (1.0 + x*(2.0 * std::log(x) - x));
}

double F3N(double x) {
   if (is_zero(x)) {
      return 8.0/105.0;
   }

   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.07)) {
      return 1.0 + d*(76/875.0 + d*(-431/2625.0 + d*(5858/42875.0
         + d*(-3561/34300.0 + d*(23/294.0 - 4381/73500.0*d)))));
   }

   return 4.0/(105.0 * pow4(d)) * (
      + (1.0 - x) * (- 97.0 * sqr(x) - 529.0 * x + 2.0)
      + 6.0 * sqr(x) * (13.0 * x + 81.0) * std::log(x)
      + 108.0 * x * (7.0 * x + 4.0) * dilog(1.0 - x)
      );
}

double F4N(double x) {
   const double PI = 3.14159265358979323846;
   const double PI2 = PI * PI;

   if (is_zero(x)) {
      return -3.0/4.0*(-9.0 + PI2);
   }

   if (is_equal(x, 1.0, 0.07)) {
      const double d = x - 1.0;

      return 1.0 + sqr(d)*(-111.0/800.0 + d*(59.0/400.0 + d*(-129.0/980.0
         + d*(177.0/1568.0 - 775.0/8064.0*d))));
   }

   return -2.25/pow3(1.0 - x) * (
      + (x + 3.0) * (x * std::log(x) + x - 1.0)
      + (6.0 * x + 2.0) * dilog(1.0 - x)
      );
}

/// Fb(1,1)
double Fb11(double x, double y) {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;

   return
      + 1.0/12.0 + (-0.05 + y1/30.)*y1
      + x1*(-0.05 + (1.0/30.0 - y1/42.0)*y1
      + x1*(1.0/30.0 + (-1.0/42.0 + y1/56.0)*y1));
}

/// Fb(x,1), x != 1, x != 0
double Fb1(double x, double y) {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;

   return (2.0 + x*(3.0 + 6.0*lx + x*(-6.0 + x)))/(6.0*x14)
      + y1*(3.0 + x*(10.0 + 12.0*lx + x*(-18.0 + x*(6.0 - x))))/(12.0*x15)
      + sqr(y1)*(12.0 + x*(65.0 + 60.0*lx + x*(-120.0 + x*(60.0 + x*(-20.0 + 3.0*x)))))/(60.*x16);
}

/// Fb(x,x), x != 1, x != 0
double Fbx(double x, double y) {
   const double x1 = x - 1.0;
   const double d = y - x;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;

   return (-5.0 - 2.0*lx + x*(4.0 - 4.0*lx + x))/(2.0*x14)
      - d*(-1.0 + x*(-9.0 - 6.0*lx + x*(9.0 - 6.0*lx + x)))/(2.0*x15*x)
      - sqr(d)*(-1.0 + x*(12.0 + x*(36.0 + 36.0*lx + x*(-44.0 + 24.0*lx - 3.0*x))))/(6.*x16*sqr(x));
}

double Fb(double x, double y) {
   if ((is_zero(x) && is_zero(y)) || is_zero(x) || is_zero(y)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.01) && is_equal(y, 1.0, 0.01)) {
      return Fb11(x, y);
   }

   if (is_equal(x, 1.0, 0.01)) {
      return Fb1(y, x);
   }

   if (is_equal(y, 1.0, 0.01)) {
      return Fb1(x, y);
   }

   if (is_equal_rel(x, y, 0.01)) {
      return Fbx(x, y);
   }

   return - (G4(x) - G4(y)) / (x - y);
}

/// Fa(1,1)
double Fa11(double x, double y) {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;

   return
      0.25 + (-0.2 + y1/6.)*y1
      + x1*(-0.2 + (1.0/6.0 - y1/7.0)*y1)
      + sqr(x1)*(1.0/6.0 + (-1.0/7.0 + y1/8.0)*y1);
}

/// Fa(x,1)
double Fa1(double x, double y) {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;

   return (-11.0 - 6.0*lx + x*(18.0 + x*(-9.0 + 2.0*x)))/(6.0*x14)
      + y1*(-25.0 - 12.0*lx + x*(48.0 + x*(-36.0 + x*(16.0 - 3.0*x))))/(12.0*x15)
      + sqr(y1)*(-137.0 - 60.0*lx + x*(300.0 + x*(-300.0 + x*(200.0 + x*(-75.0 + 12.0*x)))))/(60.0*x16);
}

/// Fa(x,x)
double Fax(double x, double y) {
   return (2. + 3.*x - 6.*sqr(x) + pow3(x) + 6.*x*std::log(x)) / (2.*pow4(-1. + x)*x)
      - (-x + y) * (-1. + 8.*x - 8.*pow3(x) + pow4(x) + 12.*sqr(x)*std::log(x))
      / (2.*std::pow(-1. + x,5)*sqr(x))
      - sqr(-x + y) * (-2. + 15.*x - 60.*sqr(x) + 20.*pow3(x) + 30.*pow4(x) - 3.*std::pow(x,5)
                       - 60.*pow3(x)*std::log(x)) / (6.*std::pow(-1. + x,6)*pow3(x));
}

double Fa(double x, double y) {
   if ((is_zero(x) && is_zero(y)) || is_zero(x) || is_zero(y)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.01) && is_equal(y, 1.0, 0.01)) {
      return Fa11(x,y);
   }

   if (is_equal(x, 1.0, 0.01)) {
      return Fa1(y,x);
   }

   if (is_equal(y, 1.0, 0.01)) {
      return Fa1(x,y);
   }

   if (is_equal_rel(x, y, 0.01)) {
      return Fax(x,y);
   }

   return - (G3(x) - G3(y)) / (x - y);
}

double G3(double x) {
   return 1.0/(2.0*pow3(x - 1.0))*((x - 1.0)*(x - 3.0) + 2.0*std::log(x));
}

double G4(double x) {
   return 1.0/(2.0*pow3(x - 1.0))*((x - 1.0)*(x + 1.0) - 2.0*x*std::log(x));
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
           - 6.*sqr(a)*sqr(c) * std::log(sqr(a)/sqr(c))
           * (6.*pow4(a) - 8.*pow3(a)*b + 3.*sqr(a)*(sqr(b) - sqr(c)) + sqr(c)*(sqr(b) + sqr(c))))
      / (6.*sqr(a)*pow4(sqr(a) - sqr(c)));
}

/// Iabc(a,a,0)
double Iaa0(double a, double b) {
   return (17.*sqr(a) - 16.*a*b + 5.*sqr(b)) / (6.*pow4(a));
}

/// Iabc(0,b,c)
double I0bc(double b, double c) {
   return std::log(sqr(b/c))/(sqr(b) - sqr(c));
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

   return ( (sqr(a * b) * std::log(sqr(a / b))
           + sqr(b * c) * std::log(sqr(b / c))
           + sqr(c * a) * std::log(sqr(c / a)))
           / ((sqr(a) - sqr(b)) * (sqr(b) - sqr(c)) * (sqr(a) - sqr(c))) );
}

/**
 * Calculates \f$f_{PS}(z)\f$, Eq (70) arXiv:hep-ph/0609168
 */
double f_PS(double z) {
   double result = 0.0;

   if (z < 0.25) {
      const double y = sqrt(1. - 4. * z);
      result = 2.0*z/y*(dilog(1.0 - 0.5*(1.0 - y)/z) - dilog(1.0 - 0.5*(1.0 + y)/z));
   } else {
      const std::complex<double> y = std::sqrt(std::complex<double>(1.0 - 4.0*z, 0.0));
      const std::complex<double> zc(z, 0.0);
      result = std::real(2.0*zc/y*(dilog(1.0 - 0.5*(1.0 - y)/zc) - dilog(1.0 - 0.5*(1.0 + y)/zc)));
   }

   return result;
}

/**
 * Calculates \f$f_S(z)\f$, Eq (71) arXiv:hep-ph/0609168
 */
double f_S(double z) {
   if (z < 0.0) {
      ERROR("f_S: z must not be negativ!");
   }

   return (2.0*z - 1.0)*f_PS(z) - 2.0*z*(2.0 + std::log(z));
}

/**
 * Calculates \f$f_{\tilde{f}}(z)\f$, Eq (72) arXiv:hep-ph/0609168
 */
double f_sferm(double z) {
   if (z < 0.0) {
      ERROR("f_sferm: z must not be negativ!");
   }

   return 0.5*z*(2.0 + std::log(z) - f_PS(z));
}

} // namespace gm2calc
