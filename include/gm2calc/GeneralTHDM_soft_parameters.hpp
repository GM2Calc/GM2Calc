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

#ifndef GeneralTHDM_SOFT_PARAMETERS_H
#define GeneralTHDM_SOFT_PARAMETERS_H

#include "GeneralTHDM_susy_parameters.hpp"

#include <iosfwd>

namespace gm2calc {

/**
 * @class GeneralTHDM_soft_parameters
 * @brief contains dimensionfull parameters of the THDM model
 *
 * Dimensionfull parameters are: scalar squared mass parameters for
 * the Higgs doublets and vacuum expectation values.
 */
class GeneralTHDM_soft_parameters : public GeneralTHDM_susy_parameters {
public:
   void print(std::ostream&) const override;

   void set_M122(double M122_) { M122 = M122_; }
   void set_M112(double M112_) { M112 = M112_; }
   void set_M222(double M222_) { M222 = M222_; }
   void set_v1(double v1_) { v1 = v1_; }
   void set_v2(double v2_) { v2 = v2_; }

   double get_M122() const { return M122; }
   double get_M112() const { return M112; }
   double get_M222() const { return M222; }
   double get_v1() const { return v1; }
   double get_v2() const { return v2; }

protected:
   double M122{0.0};
   double M112{0.0};
   double M222{0.0};
   double v1{0.0};
   double v2{0.0};
};

std::ostream& operator<<(std::ostream&, const GeneralTHDM_soft_parameters&);

} // namespace gm2calc

#endif
