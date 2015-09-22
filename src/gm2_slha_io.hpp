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

#ifndef GM2_SLHA_IO_H
#define GM2_SLHA_IO_H

#include <string>
#include <iosfwd>
#include <functional>
#include <Eigen/Core>
#include <boost/format.hpp>
#include "slhaea.h"
#include "gm2_error.hpp"
#include "numerics2.hpp"

namespace gm2calc {

class MSSMNoFV_onshell;
struct MSSMNoFV_onshell_physical;

   namespace {
      /// SLHA number formatter
      const boost::format number_formatter("         %16.8E   # %s\n");
      /// SLHA scale formatter
      const boost::format scale_formatter("%9.8E");
      /// SLHA line formatter for the one-element entries (HMIX, GAUGE, MSOFT, ...)
      const boost::format single_element_formatter(" %5d   %16.8E   # %s\n");
      /// SLHA line formatter for the SPINFO block entries
      const boost::format spinfo_formatter(" %5d   %s\n");
   }

#define FORMAT_ELEMENT(pdg,value,name)                                  \
   boost::format(single_element_formatter) % pdg % value % name
#define FORMAT_SCALE(n)                                                 \
   boost::format(scale_formatter) % n
#define FORMAT_NUMBER(n,str)                                            \
   boost::format(number_formatter) % n % str
#define FORMAT_SPINFO(n,str)                                            \
   boost::format(spinfo_formatter) % n % str

/**
 * @class Config_options
 * @brief configuration for the calculation of \f$a_\mu\f$
 */
struct Config_options {
   enum E_output_format : unsigned {
      Minimal = 0, Detailed = 1, NMSSMTools = 2, SPheno = 3 };

   E_output_format output_format = Minimal; ///< output format
   unsigned loop_order = 2;      ///< loop order
   bool tanb_resummation = true; ///< tan(beta) resummation
   bool force_output = false;    ///< print output even if error occured
   bool verbose_output = false;  ///< print additional information
};

/**
 * @class GM2_slha_io
 * @brief class for reading input files and writing SLHA output files
 */
class GM2_slha_io {
public:
   typedef std::function<void(int, double)> Tuple_processor;
   enum Position { front, back };

   void clear();

   // reading functions
   const SLHAea::Coll& get_data() const { return data; }
   bool block_exists(const std::string&) const;
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   void read_block(const std::string&, const Tuple_processor&, double scale = 0) const;
   template <class Derived>
   void read_block(const std::string&, Eigen::MatrixBase<Derived>&, double scale = 0) const;
   double read_entry(const std::string&, int, double scale = 0) const;
   double read_scale(const std::string&) const;

   // writing functions
   void set_data(const SLHAea::Coll& data_) { data = data_; }
   void set_block(const std::ostringstream&, Position position = back);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& = std::cout);
   void fill_block_entry(const std::string&, unsigned, double, const std::string&);
   void fill_block_entry(const std::string&, unsigned, const std::string&);

private:
   SLHAea::Coll data;          ///< SHLA data
   template <class Scalar>
   static Scalar convert_to(const std::string&); ///< convert string
   static std::string to_lower(const std::string&); ///< string to lower case
   static bool at_scale(const SLHAea::Block&, double); ///< check block scale
};

template <class Scalar>
Scalar GM2_slha_io::convert_to(const std::string& str)
{
   Scalar value;
   try {
      value = SLHAea::to<Scalar>(str);
   }  catch (const boost::bad_lexical_cast& error) {
      const std::string msg("cannot convert string \"" + str + "\" to "
                            + typeid(Scalar).name());
      throw EReadError(msg);
   }
   return value;
}

/**
 * Fills a matrix from a SLHA block
 *
 * @param block_name block name
 * @param matrix matrix to be filled
 * @param scale (or 0 if scale should be ignored)
 */
template <class Derived>
void GM2_slha_io::read_block(const std::string& block_name,
                             Eigen::MatrixBase<Derived>& matrix,
                             double scale) const
{
   SLHAea::Coll::const_iterator block =
      data.find(data.cbegin(), data.cend(), block_name);

   const int cols = matrix.cols(), rows = matrix.rows();

   while (block != data.cend()) {
      if (!at_scale(*block, scale))
         continue;

      for (SLHAea::Block::const_iterator line = block->cbegin(),
              end = block->cend(); line != end; ++line) {
         if (!line->is_data_line())
            continue;

         if (cols == 1) {
            // vector
            if (line->size() >= 2) {
               const int i = convert_to<int>((*line)[0]) - 1;
               if (0 <= i && i < rows) {
                  const double value = convert_to<double>((*line)[1]);
                  matrix(i,0) = value;
               }
            }
         } else {
            // martix
            if (line->size() >= 3) {
               const int i = convert_to<int>((*line)[0]) - 1;
               const int k = convert_to<int>((*line)[1]) - 1;
               if (0 <= i && i < rows && 0 <= k && k < cols) {
                  const double value = convert_to<double>((*line)[2]);
                  matrix(i,k) = value;
               }
            }
         }
      }

      ++block;
      block = data.find(block, data.cend(), block_name);
   }
}

/// read model parameters (GM2Calc input format)
void fill_gm2calc(const GM2_slha_io&, MSSMNoFV_onshell&);

/// read model parameters (SLHA input format)
void fill_slha(const GM2_slha_io&, MSSMNoFV_onshell&);

/// read configuration
void fill(const GM2_slha_io&, Config_options&);

} // namespace gm2calc

#endif
