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

#ifndef GM2_FFUNCTIONS_H
#define GM2_FFUNCTIONS_H

namespace gm2calc {

/// \f$F_1^C(x)\f$, Eq (54) arXiv:hep-ph/0609168
double F1C(double);
/// \f$F_2^C(x)\f$, Eq (55) arXiv:hep-ph/0609168
double F2C(double);
double F3C(double);
double F4C(double);
/// \f$F_1^N(x)\f$, Eq (52) arXiv:hep-ph/0609168
double F1N(double);
/// \f$F_2^N(x)\f$, Eq (53) arXiv:hep-ph/0609168
double F2N(double);
double F3N(double);
double F4N(double);
/// \f$F_a(x)\f$, Eq (6.3a) arXiv:1311.1775
double Fa(double, double);
/// \f$F_b(x)\f$, Eq (6.3b) arXiv:1311.1775
double Fb(double, double);
/// \f$G_3(x)\f$, Eq (6.4a) arXiv:1311.1775
double G3(double);
/// \f$G_4(x)\f$, Eq (6.4b) arXiv:1311.1775
double G4(double);
/// \f$H_2(x,y)\f$, Eq (34) of arxiv:0901.2065
double H2(double, double);
double Iabc(double, double, double);
/// \f$f_{PS}(z)\f$, Eq (70) arXiv:hep-ph/0609168
double f_PS(double);
/// \f$f_S(z)\f$, Eq (71) arXiv:hep-ph/0609168
double f_S(double);
/// \f$f_{\tilde{f}}(z)\f$, Eq (72) arXiv:hep-ph/0609168
double f_sferm(double);

double cube(double);
double quad(double);
int sign(double);
double signed_abs_sqrt(double);
double signed_sqr(double);
template <typename T> T sqr(T x) { return x*x; }

} // namespace gm2calc

#endif
