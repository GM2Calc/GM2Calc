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

#ifndef GM2_ERROR_H
#define GM2_ERROR_H

#include "error.hpp"

namespace flexiblesusy {
namespace gm2calc {

class EInvalidInput : public Error {
public:
   explicit EInvalidInput(const std::string& message_) : message(message_) {}
   virtual ~EInvalidInput() {}
   virtual std::string what() const { return message; }
private:
   std::string message;
};

} // namespace gm2calc
} // namespace flexiblesusy

#endif
