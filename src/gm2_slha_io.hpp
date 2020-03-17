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

#include "gm2calc/gm2_error.hpp"

#include "slhaea.h"

#include <cmath>
#include <string>
#include <iostream>

#include <Eigen/Core>
#include <boost/format.hpp>

namespace gm2calc {

struct Config_options;
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
 * @class GM2_slha_io
 * @brief class for reading input files and writing SLHA output files
 */
class GM2_slha_io {
public:
   typedef std::function<void(int, double)> Tuple_processor;

   void clear();

   // reading functions
   const SLHAea::Coll& get_data() const { return data; }
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   void read_block(const std::string&, const Tuple_processor&, double scale = 0) const;
   template <class Derived>
   void read_block(const std::string&, Eigen::MatrixBase<Derived>&, double scale = 0) const;
   double read_scale(const std::string&) const;

   // writing functions
   void set_data(const SLHAea::Coll& data_) { data = data_; }
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream&);
   void fill_block_entry(const std::string&, unsigned, double, const std::string&);
   void fill_block_entry(const std::string&, unsigned, const std::string&);

   /// read model parameters (GM2Calc input format)
   void fill_gm2calc(MSSMNoFV_onshell&) const;

   /// read model parameters (SLHA input format)
   void fill_slha(MSSMNoFV_onshell&) const;

   /// read configuration
   void fill(Config_options&) const;

private:
   SLHAea::Coll data;          ///< SHLA data
   template <class Scalar>
   static Scalar convert_to(const std::string&); ///< convert string
   static bool is_at_scale(const SLHAea::Block&, double, double eps = 0.01); ///< check block scale
   static double read_scale(const SLHAea::Block&); ///< read scale from block
   static void read_block(const SLHAea::Block&, const Tuple_processor&); ///< read block
   template <class Derived>
   static void read_block(const SLHAea::Block&, Eigen::MatrixBase<Derived>&);

   void fill_scale(MSSMNoFV_onshell&) const;
   void fill_alpha_from_gm2calcinput(MSSMNoFV_onshell&) const;
   void fill_from_A(MSSMNoFV_onshell&) const;
   void fill_from_gm2calcinput(MSSMNoFV_onshell&) const;
   void fill_from_hmix(MSSMNoFV_onshell&) const;
   void fill_from_mass(MSSMNoFV_onshell_physical&) const;
   void fill_from_msoft(MSSMNoFV_onshell&) const;
   void fill_from_sminputs(MSSMNoFV_onshell&) const;
};

template <class Scalar>
Scalar GM2_slha_io::convert_to(const std::string& str)
{
   Scalar value;
   try {
      value = SLHAea::to<Scalar>(str);
      if (!std::isfinite(value)) {
         throw boost::bad_lexical_cast();
      }
   }  catch (const boost::bad_lexical_cast& error) {
      throw EReadError("non-numeric input");
   }
   return value;
}

/**
 * Fills a matrix from an SLHA block
 *
 * @param block the block
 * @param matrix matrix to be filled
 */
template <class Derived>
void GM2_slha_io::read_block(const SLHAea::Block& block, Eigen::MatrixBase<Derived>& matrix)
{
   const int cols = matrix.cols(), rows = matrix.rows();

   // vector
   if (cols == 1) {
      for (const auto& line : block) {
         if (line.is_data_line() && line.size() >= 2) {
            const int i = convert_to<int>(line[0]) - 1;
            if (0 <= i && i < rows) {
               const double value = convert_to<double>(line[1]);
               matrix(i) = value;
            }
         }
      }
   }

   // matrix
   if (cols > 1) {
      for (const auto& line: block) {
         if (line.is_data_line() && line.size() >= 3) {
            const int i = convert_to<int>(line[0]) - 1;
            const int k = convert_to<int>(line[1]) - 1;
            if (0 <= i && i < rows && 0 <= k && k < cols) {
               const double value = convert_to<double>(line[2]);
               matrix(i,k) = value;
            }
         }
      }
   }
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

   while (block != data.cend()) {
      if (is_at_scale(*block, scale)) {
         GM2_slha_io::read_block(*block, matrix);
      }

      ++block;
      block = data.find(block, data.cend(), block_name);
   }
}

} // namespace gm2calc

#endif
