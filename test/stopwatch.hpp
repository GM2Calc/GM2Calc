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

#ifndef STOPWATCH_HPP
#define STOPWATCH_HPP

#include <chrono>

namespace gm2calc {

class Stopwatch {
private:
   using microseconds_t = std::chrono::duration<int,std::micro>;

   microseconds_t::rep get_ticks() const {
      microseconds_t duration(std::chrono::duration_cast<microseconds_t>(
                                 stop_point - start_point));
      return duration.count();
   }

public:
   void start() {
      start_point = std::chrono::high_resolution_clock::now();
   }

   void stop() {
      stop_point = std::chrono::high_resolution_clock::now();
   }

   double get_time_in_seconds() const {
      return get_ticks() * 0.000001;
   }

   double get_time_in_milliseconds() const {
      return get_ticks() * 0.001;
   }

private:
   std::chrono::high_resolution_clock::time_point start_point, stop_point;
};

} // namespace gm2calc

#endif
