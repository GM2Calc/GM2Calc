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

#include <stdexcept>

namespace gm2calc {

class Error : public std::runtime_error {
public:
   Error(const char* msg) : std::runtime_error(msg) {}
};

/**
 * @class ESetupError
 * @brief Spectrum generator was not setup correctly
 */
class ESetupError : public Error {
public:
   ESetupError(const char* msg) : Error(msg) {}
};

class EInvalidInput : public Error {
public:
   EInvalidInput(const char* msg) : Error(msg) {}
};

class EPhysicalProblem : public Error {
public:
   EPhysicalProblem(const char* msg) : Error(msg) {}
};

class EReadError : public Error {
public:
   EReadError(const char* msg) : Error(msg) {}
};

} // namespace gm2calc

#endif