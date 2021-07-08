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

#ifndef GM2_GeneralTHDM_HPP
#define GM2_GeneralTHDM_HPP

#include "GeneralTHDM_mass_eigenstates.hpp"

#include <iosfwd>

#include <Eigen/Core>

namespace gm2calc {

/**
 * @class GeneralTHDM
 * @brief Contains routines to determine the THDM parameters
 */
class GeneralTHDM : public GeneralTHDM_mass_eigenstates {
public:
   enum class Yukawa_scheme {
      type_1, type_2, type_X, type_Y
   };

   virtual ~GeneralTHDM() = default;

   double get_zeta_u() const;
   double get_zeta_d() const;
   double get_zeta_l() const;

   void set_yukawa_scheme(Yukawa_scheme ys) { yukawa_scheme = ys; }

private:
   Yukawa_scheme yukawa_scheme{Yukawa_scheme::type_2};
};

/// streaming operator
std::ostream& operator<<(std::ostream&, const GeneralTHDM&);

} // namespace gm2calc

#endif
