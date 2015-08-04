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

namespace flexiblesusy {

namespace {

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

} // anonymous namespace

namespace gm2calc {

double sqr(double x) { return x*x; }
double cube(double x) { return x*x*x; }
double quad(double x) { return sqr(x)*sqr(x); }
int sign(double x) { return x < 0 ? -1 : 1; }
double signed_sqr(double x) { return sign(x) * x * x; }
double signed_abs_sqrt(double x) {
   return sign(x) * std::sqrt(std::abs(x));
}

double F1C(double x) {
   if(is_equal(x, 1.))
      ERROR("F1C: x must not be 1 !");

   return 2. / quad(1. - x) * (2. + 3. * x - 6. * sqr(x)
                           + cube(x) + 6. * x * log(x));
}

double F2C(double x) {
   if(is_equal(x, 1.))
      ERROR("F2C: x must not be 1 !");

   return 3. / (2. * cube(1. - x)) * (- 3. + 4. * x - sqr(x) - 2. * log(x));
}

double F3C(double x) {
   if(is_equal(x, 1.))
      ERROR("F3C: x must not be 1 !");

   return ( 4. / (141. * quad(1. - x)) * ((1. - x) * (151. * sqr(x) - 335. * x + 592.)
             + 6. * (21. * cube(x) - 108. * sqr(x) - 93. * x + 50.) * log(x)
             - 54. * x * (sqr(x) - 2. * x - 2.) * sqr(log(x))
             - 108. * x * (sqr(x) - 2. * x + 12.) * dilog(1.- x)) );
}

double F4C(double x) {
   if(is_equal(x, 1.))
      ERROR("F4C: x must not be 1 !");

   return ( - 9. / (122. * cube(1. - x)) * (8. * (sqr(x) - 3. * x + 2.)
             +  (11. * sqr(x) - 40. * x + 5.) * log(x)
             - 2. * (sqr(x) - 2. * x - 2.) * sqr(log(x))
             - 4. * (sqr(x) - 2. * x + 9.) * dilog(1.- x)) );
}

double F1N(double x) {
   if(is_equal(x, 1.))
      ERROR("F1N: x must not be 1 !");

   return 2. / quad(1. - x) * (1. - 6. * x + 3. * sqr(x)
                          + 2. * cube(x) - 6. * sqr(x) * log(x));
}

double F2N(double x) {
   if(is_equal(x, 1.))
      ERROR("F2N: x must not be 1 !");

   return 3. / cube(1. - x) * (1. - sqr(x) + 2. * x * log(x));
}

double F3N(double x) {
   if(is_equal(x, 1.))
      ERROR("F3N: x must not be 1 !");

   return 4. / (105. * quad(1. - x)) * ((1. - x) * (- 97. * sqr(x) - 529. * x + 2.)
            + 6. * sqr(x) * (13. * x + 81.) * log(x)
            + 108. * x * (7. * x + 4.) * dilog(1. - x));
}

double F4N(double x) {
   if(is_equal(x, 1.))
      ERROR("F4N: x must not be 1 !");

   return - 2.25 / cube(1. - x) * ((x + 3.) * (x * log(x) + x - 1.)
                                  + (6. * x + 2.) * dilog(1. - x));
}

double Fa(double x, double y) {
   if(is_equal(x, y))
      ERROR("Fa: x must not be equal y!");

   return - (G3(x) - G3(y)) / (x - y);
}

double Fb(double x, double y) {
   if(is_equal(x, y))
      ERROR("Fb: x must not be equal y!");

   return - (G4(x) - G4(y)) / (x - y);
}

double G3(double x) {
   if(is_equal(x, 1.))
      ERROR("G3: x must not be 1 !");

   return 1. / (2. * cube(x - 1.)) * ((x - 1.) * (x - 3.) + 2. * log(x));
}

double G4(double x) {
   if(is_equal(x, 1.))
      ERROR("G4: x must not be 1 !");

   return 1. / (2. * cube(x - 1.)) * ((x - 1.) * (x + 1.) - 2. * x * log(x));
}

/**
 * \f$H_2(x,y)\f$ function from Eq. (34) of arxiv:0901.2065 .
 *
 * The \f$H_2(x,y)\f$ function is related to \f$I(a,b,c)\f$ by
 * \f$\frac{1}{z} H_2(\frac{x}{z},\frac{y}{z}) = - I(x,y,z)\f$.
 */
double H2(double x, double y) {
   if (is_equal(x,1.) || is_equal(y,1.) || is_equal(x,y))
      ERROR("H2(" << x << "," << y << ") is not well-defined!");

   return x * log(x) / ((1-x)*(x-y)) + y * log(y) / ((1-y)*(y-x));
}

double Iaac(double a, double c) {
   const double va2 = sqr(a/c);
   return (va2 - 1. - log(va2))/(sqr(c) * sqr(va2 - 1.));
}

double Iabb(double a, double b) {
   const double va = a/b;
   const double va2 = sqr(va);
   return va*(1. - va2 + va2*log(va2))/(a*b*sqr(va2 - 1.));
}

double Iabc(double a, double b, double c) {
   if (is_equal(a,b) && is_equal(b,c))
      return 0.5 / sqr(a);

   if (is_equal(a,b))
      return Iaac(a,c);

   if (is_equal(b,c))
      return Iabb(a,b);

   if (is_equal(a,c))
      return Iabb(b,a);

   return ( (sqr(a * b) * log(sqr(a / b))
           + sqr(b * c) * log(sqr(b / c))
           + sqr(c * a) * log(sqr(c / a)))
           / ((sqr(a) - sqr(b)) * (sqr(b) - sqr(c)) * (sqr(a) - sqr(c))) );
}

/**
 * Calculates \f$f_{PS}(z)\f$, Eq. (70) arXiv:hep-ph/0609168
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
 * Calculates \f$f_S(z)\f$, Eq. (71) arXiv:hep-ph/0609168
 */
double f_S(double z) {
   if(z < 0.)
      ERROR("f_S: z must not be negativ!");

   return (2. * z - 1.) * f_PS(z) - 2. * z * (2. + log(z));
}

/**
 * Calculates \f$f_{\tilde{f}}(z)\f$, Eq. (72) arXiv:hep-ph/0609168
 */
double f_sferm(double z) {
   if(z < 0.)
      ERROR("f_sferm: z must not be negativ!");

   return 0.5 * z * (2. + log(z) - f_PS(z));
}

} // namespace gm2calc
} // namespace flexiblesusy
