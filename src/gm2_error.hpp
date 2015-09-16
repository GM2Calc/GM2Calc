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

namespace gm2calc {

class Error {
public:
   virtual ~Error() {}
   virtual std::string what() const = 0;
};

/**
 * @class ESetupError
 * @brief Spectrum generator was not setup correctly
 */
class ESetupError : public Error {
public:
   explicit ESetupError(const std::string& message_) : message(message_) {}
   virtual ~ESetupError() {}
   virtual std::string what() const { return message; }
private:
   std::string message;
};

class EInvalidInput : public Error {
public:
   explicit EInvalidInput(const std::string& message_) : message(message_) {}
   virtual ~EInvalidInput() {}
   virtual std::string what() const { return message; }
private:
   std::string message;
};

class EPhysicalProblem : public Error {
public:
   explicit EPhysicalProblem(const std::string& message_) : message(message_) {}
   virtual ~EPhysicalProblem() {}
   virtual std::string what() const { return message; }
private:
   std::string message;
};

class EReadError : public Error {
public:
   EReadError(const std::string& message_) : message(message_) {}
   virtual ~EReadError() {}
   virtual std::string what() const { return message; }
private:
   std::string message;
};

} // namespace gm2calc

#endif
