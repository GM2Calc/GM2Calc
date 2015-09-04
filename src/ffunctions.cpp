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
#include "dilog.h"
#include "numerics2.hpp"

#include <iostream>
#include <cmath>
#include <complex>

#define ERROR(message) std::cerr << "Error: " << message << '\n';

namespace gm2calc {

using namespace flexiblesusy;

double dilog(double x) {
   // The DiLogarithm function
   // Code translated by R.Brun from CERNLIB DILOG function C332

   const double PI = M_PI;
   const double HF  = 0.5;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double PI12 = PI2/12;
   const double C[20] = {0.42996693560813697, 0.40975987533077105,
     -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
      0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
     -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
      0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
     -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
      0.00000000000000093,-0.00000000000000014, 0.00000000000000002};

   double T,H,Y,S,A,ALFA,B1,B2,B0;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= log(-T);
           B2= log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = log(-T);
           A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5) {
           Y = -(1+T)/T;
           S = 1;
           A = log(-T);
           A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i=19;i>=0;i--){
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}

std::complex<double> dilog(const std::complex<double>& x) {
   double a = x.real(), b = x.imag();
   double ansreal = 0., ansimag = 0.;

   dilogwrap_(&a, &b, &ansreal, &ansimag);

   return std::complex<double>(ansreal, ansimag);
}

double sqr(double x) { return x*x; }
double cube(double x) { return x*x*x; }
double quad(double x) { return sqr(x)*sqr(x); }
int sign(double x) { return x < 0 ? -1 : 1; }
double signed_sqr(double x) { return sign(x) * x * x; }
double signed_abs_sqrt(double x) {
   return sign(x) * std::sqrt(std::abs(x));
}

double F1C(double x) {
   if (is_equal(x, 1., 0.05)) {
      return 1. - 0.6*(-1 + x) + 0.4*sqr(-1 + x) - 2/7.*cube(-1 + x)
         + 3/14.*quad(-1 + x) - 1/6.*std::pow(-1 + x,5)
         + 2/15.*std::pow(-1 + x,6);
   }

   return 2. / quad(1. - x) * (2. + 3. * x - 6. * sqr(x)
                               + cube(x) + 6. * x * log(x));
}

double F2C(double x) {
   if (is_equal(x, 1., 0.03)) {
      return 1. - 0.75*(-1 + x) + 0.6*sqr(-1 + x) - 0.5*cube(-1 + x)
         + 3/7.*quad(-1 + x) - 0.375*std::pow(-1 + x,5)
         + 1/3.*std::pow(-1 + x,6);
   }

   return 3. / (2. * cube(1. - x)) * (- 3. + 4. * x - sqr(x) - 2. * log(x));
}

double F3C(double x) {
   if (is_equal(x, 1., 0.07)) {
      return 1. + 1059/1175.*(-1 + x)
         - 4313/3525.*sqr(-1 + x)
         + 70701/57575.*cube(-1 + x)
         - 265541/230300.*quad(-1 + x)
         + 48919/46060.*std::pow(-1 + x,5)
         - 80755/82908.*std::pow(-1 + x,6);
   }

   return ( 4. / (141. * quad(1. - x)) * ((1. - x) * (151. * sqr(x) - 335. * x + 592.)
             + 6. * (21. * cube(x) - 108. * sqr(x) - 93. * x + 50.) * log(x)
             - 54. * x * (sqr(x) - 2. * x - 2.) * sqr(log(x))
             - 108. * x * (sqr(x) - 2. * x + 12.) * dilog(1.- x)) );
}

double F4C(double x) {
   if (is_equal(x, 1., 0.07)) {
      return 1. - 45/122.*(-1 + x)
         + 941/6100.*sqr(-1 + x)
         - 17/305.*cube(-1 + x)
         + 282/74725.*quad(-1 + x)
         + 177/6832.*std::pow(-1 + x,5)
         - 47021/1076040.*std::pow(-1 + x,6);
   }

   return ( - 9. / (122. * cube(1. - x)) * (8. * (sqr(x) - 3. * x + 2.)
             +  (11. * sqr(x) - 40. * x + 5.) * log(x)
             - 2. * (sqr(x) - 2. * x - 2.) * sqr(log(x))
             - 4. * (sqr(x) - 2. * x + 9.) * dilog(1.- x)) );
}

double F1N(double x) {
   if (is_equal(x, 1., 0.06)) {
      return 1. - 0.4*(-1 + x) + 0.2*sqr(-1 + x) - 4/35.*cube(-1 + x)
         + 1/14.*quad(-1 + x) - 1/21.*std::pow(-1 + x,5)
         + 1/30.*std::pow(-1 + x,6);
   }

   return 2. / quad(1. - x) * (1. - 6. * x + 3. * sqr(x)
                               + 2. * cube(x) - 6. * sqr(x) * log(x));
}

double F2N(double x) {
   if (is_equal(x, 1., 0.04)) {
      return 1. - 0.5*(-1 + x) + 0.3*sqr(-1 + x) - 0.2*cube(-1 + x)
         + 1/7.*quad(-1 + x) - 3/28.*std::pow(-1 + x,5)
         + 1/12.*std::pow(-1 + x,6);
   }

   return 3. / cube(1. - x) * (1. - sqr(x) + 2. * x * log(x));
}

double F3N(double x) {
   if (is_equal(x, 1., 0.07)) {
      return 1. + 76/875.*(-1 + x)
         - 431/2625.*sqr(-1 + x)
         + 5858/42875.*cube(-1 + x)
         - 3561/34300.*quad(-1 + x)
         + 23/294.*std::pow(-1 + x,5)
         - 4381/73500.*std::pow(-1 + x,6);
   }

   return 4. / (105. * quad(1. - x)) * ((1. - x) * (- 97. * sqr(x) - 529. * x + 2.)
            + 6. * sqr(x) * (13. * x + 81.) * log(x)
            + 108. * x * (7. * x + 4.) * dilog(1. - x));
}

double F4N(double x) {
   if (is_equal(x, 1., 0.07)) {
      return 1. - 0.13875*sqr(-1 + x)
         + 0.1475*cube(-1 + x)
         - 129/980.*quad(-1 + x)
         + 177/1568.*std::pow(-1 + x,5)
         - 775/8064.*std::pow(-1 + x,6);
   }

   return - 2.25 / cube(1. - x) * ((x + 3.) * (x * log(x) + x - 1.)
                                   + (6. * x + 2.) * dilog(1. - x));
}

/// Fb(1,1)
double Fb11(double x, double y) {
   return (293. - 216.*y + 63.*sqr(y) + sqr(x)*(63. - 50.*y + 15.*sqr(y))
           - 2.*x*(108. - 84.*y + 25.*sqr(y))) / 840.;
}

/// Fb(x,1)
double Fb1(double x, double y) {
   return (2. + 3.*x - 6.*sqr(x) + cube(x) + 6.*x*log(x)) / (6.*quad(-1. + x))
      + (-1. + y) * (3. + 10.*x - 18.*sqr(x) + 6.*cube(x) - quad(x) + 12.*x*log(x))
      / (12.*std::pow(-1. + x,5))
      + sqr(-1. + y) * (12. + 65.*x - 120.*sqr(x) + 60.*cube(x) - 20.*quad(x)
                        + 3.*std::pow(x,5) + 60.*x*log(x)) / (60.*std::pow(-1. + x,6));
}

/// Fb(x,x)
double Fbx(double x, double y) {
   return (-5. + 4.*x + sqr(x) - 2.*log(x) - 4.*x*log(x)) / (2.*quad(-1. + x))
      - (-x + y) * (-1. - 9.*x + 9.*sqr(x) + cube(x) - 6.*x*log(x) - 6.*sqr(x)*log(x))
      / (2.*std::pow(-1. + x,5)*x)
      - sqr(-x + y) * (-1. + 12.*x + 36.*sqr(x) - 44.*cube(x) - 3.*quad(x) + 36.*sqr(x)*log(x)
                       + 24.*cube(x)*log(x)) / (6.*std::pow(-1. + x,6)*sqr(x));
}

double Fb(double x, double y) {
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
   return (-11. + 18.*x - 9.*sqr(x) + 2.*cube(x) - 6.*log(x)) / (6.*quad(-1. + x))
      + (-1. + y) * (-25. + 48.*x - 36.*sqr(x) + 16.*cube(x) - 3.*quad(x) - 12.*log(x))
      / (12.*std::pow(-1. + x,5))
      + sqr(-1. + y) * (-137. + 300.*x - 300.*sqr(x) + 200.*cube(x) - 75.*quad(x)
                        + 12.*std::pow(x,5) - 60.*log(x)) / (60.*std::pow(-1. + x,6));
}

/// Fa(x,x)
double Fax(double x, double y) {
   return (2. + 3.*x - 6.*sqr(x) + cube(x) + 6.*x*log(x)) / (2.*quad(-1. + x)*x)
      - (-x + y) * (-1. + 8.*x - 8.*cube(x) + quad(x) + 12.*sqr(x)*log(x))
      / (2.*std::pow(-1. + x,5)*sqr(x))
      - sqr(-x + y) * (-2. + 15.*x - 60.*sqr(x) + 20.*cube(x) + 30.*quad(x) - 3.*std::pow(x,5)
                       - 60.*cube(x)*log(x)) / (6.*std::pow(-1. + x,6)*cube(x));
}

double Fa(double x, double y) {
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
   return 1. / (2. * cube(x - 1.)) * ((x - 1.) * (x - 3.) + 2. * log(x));
}

double G4(double x) {
   return 1. / (2. * cube(x - 1.)) * ((x - 1.) * (x + 1.) - 2. * x * log(x));
}

/**
 * \f$H_2(x,y)\f$ function from Eq (34) of arxiv:0901.2065 .
 *
 * The \f$H_2(x,y)\f$ function is related to \f$I(a,b,c)\f$ by
 * \f$\frac{1}{z^2} H_2(\frac{x^2}{z^2},\frac{y^2}{z^2}) = - I(x,y,z)\f$.
 */
double H2(double x, double y) {
   if (is_equal(x, 1., 1e-8) && is_equal(y, 1., 1e-8))
      return -0.5;

   if (is_equal(x, 1., 1e-8))
      return (-1. + y - y*log(y))/sqr(-1. + y);

   if (is_equal(y, 1., 1e-8))
      return (-1. + x - x*log(x))/sqr(-1. + x);

   if (is_equal(x, y, 1e-8))
      return (1 - y + log(y))/sqr(-1 + y);

   return x * log(x) / ((1-x)*(x-y)) + y * log(y) / ((1-y)*(y-x));
}

/// Iabc(a,a,a)
double Iaaa(double a, double b, double c) {
   return (151.*quad(a) + 13.*sqr(b)*sqr(c) - 128.*cube(a)*(b + c) - 40.*a*b*c*(b + c)
           + sqr(a)*(37.*sqr(b) + 128.*b*c + 37.*sqr(c))) / (60.*std::pow(a,6));
}

/// Iabc(a,a,c)
double Iaac(double a, double b, double c) {
   return ((sqr(a) - sqr(c))
           * (17.*std::pow(a,6) - 16.*std::pow(a,5)*b - 40.*cube(a)*b*sqr(c)
              + 8.*a*b*quad(c) - sqr(b)*quad(c) + quad(a)*(5.*sqr(b) + 8.*sqr(c))
              + sqr(a)*(20.*sqr(b)*sqr(c) - quad(c)))
           - 6.*sqr(a)*sqr(c) * log(sqr(a)/sqr(c))
           * (6.*quad(a) - 8.*cube(a)*b + 3.*sqr(a)*(sqr(b) - sqr(c)) + sqr(c)*(sqr(b) + sqr(c))))
      / (6.*sqr(a)*quad(sqr(a) - sqr(c)));
}

double Iabc(double a, double b, double c) {
   if (is_equal_rel(std::abs(a), std::abs(b), 0.01) && is_equal_rel(std::abs(a), std::abs(c), 0.01))
      return Iaaa(std::abs(a),std::abs(b),std::abs(c));

   if (is_equal_rel(std::abs(a), std::abs(b), 0.01))
      return Iaac(std::abs(a),std::abs(b),c);

   if (is_equal_rel(std::abs(b), std::abs(c), 0.01))
      return Iaac(std::abs(b),std::abs(c),a);

   if (is_equal_rel(std::abs(a), std::abs(c), 0.01))
      return Iaac(std::abs(a),std::abs(c),b);

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
   if(z < 0.25) {
      double y = sqrt(1. - 4. * z);
      result = 2. * z / y * (dilog(1. - 0.5 * (1. - y) / z) - dilog(1. - 0.5 * (1. + y) / z));
   } else {
      std::complex<double> y = sqrt(std::complex<double>(1. - 4. * z, 0.));
      std::complex<double> zc(z, 0.);
      result = real(2. * zc / y * (dilog(1. - 0.5 * (1. - y) / zc) - dilog(1. - 0.5 * (1. + y) / zc)));
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
