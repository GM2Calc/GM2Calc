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

#include "gm2_slha_io.hpp"
#include "ffunctions.hpp"
#include "MSSMNoFV_onshell.hpp"

#include <cmath>
#include <fstream>
#include <limits>

#include <Eigen/Core>

namespace gm2calc {

using namespace flexiblesusy;

#define ERROR(message) std::cerr << "Error: " << message << '\n';
#define WARNING(message) std::cerr << "Warning: " << message << '\n';

namespace {

   void process_gm2calcconfig_tuple(Config_options&, int, double);
   void process_gm2calcinput_tuple(MSSMNoFV_onshell&, int, double);
   void process_fermion_sminputs_tuple(MSSMNoFV_onshell_physical&, int, double);
   void process_mass_tuple(MSSMNoFV_onshell_physical&, int, double);
   void process_msoft_tuple(MSSMNoFV_onshell&, int, double);

} // anonymous namespace

/**
 * @brief reads from source
 *
 * If source is "-", then read_from_stream() is called.  Otherwise,
 * read_from_file() is called.
 *
 * @param source string that specifies the source
 */
void GM2_slha_io::read_from_source(const std::string& source)
{
   if (source == "-")
      read_from_stream(std::cin);
   else
      read_from_file(source);
}

/**
 * @brief opens SLHA input file and reads the content
 * @param file_name SLHA input file name
 */
void GM2_slha_io::read_from_file(const std::string& file_name)
{
   std::ifstream ifs(file_name);
   if (ifs.good()) {
      data.clear();
      data.read(ifs);
   } else {
      std::ostringstream msg;
      msg << "cannot read SLHA file: \"" << file_name << "\"";
      throw EReadError(msg.str());
   }
}

/**
 * @brief reads SLHA data from a stream
 * @param istr input stream
 */
void GM2_slha_io::read_from_stream(std::istream& istr)
{
   data.read(istr);
}

double GM2_slha_io::read_entry(const std::string& block_name, int key,
                               double scale) const
{
   SLHAea::Coll::const_iterator block =
      data.find(data.cbegin(), data.cend(), block_name);

   double entry = 0.;
   const SLHAea::Block::key_type keys(1, std::to_string(key));
   SLHAea::Block::const_iterator line;

   while (block != data.cend()) {
      if (at_scale(*block, scale)) {
         line = block->find(keys);

         if (line != block->end() && line->is_data_line() && line->size() > 1)
            entry = convert_to<double>(line->at(1));
      }

      ++block;
      block = data.find(block, data.cend(), block_name);
   }

   return entry;
}

/**
 * Reads scale definition from SLHA block.
 *
 * @param block_name block name
 *
 * @return scale (or 0 if no scale is defined)
 */
double GM2_slha_io::read_scale(const std::string& block_name) const
{
   if (!block_exists(block_name))
      return 0.;

   double scale = 0.;

   for (SLHAea::Block::const_iterator line = data.at(block_name).cbegin(),
        end = data.at(block_name).cend(); line != end; ++line) {
      if (!line->is_data_line()) {
         if (line->size() > 3 &&
             to_lower((*line)[0]) == "block" && (*line)[2] == "Q=")
            scale = convert_to<double>((*line)[3]);
         break;
      }
   }

   return scale;
}

bool GM2_slha_io::block_exists(const std::string& block_name) const
{
   return data.find(block_name) != data.cend();
}

/**
 * Returns true if the block scale after Q= matches \a scale, false
 * otherwise.  If scale == 0, the functions returns true.
 */
bool GM2_slha_io::at_scale(const SLHAea::Block& block, double scale)
{
   if (flexiblesusy::is_zero(scale))
      return true;

   for (SLHAea::Block::const_iterator line = block.cbegin(),
           end = block.cend(); line != end; ++line) {
      // check scale from block definition matches argument
      if (!line->is_data_line() && line->size() > 3 &&
          to_lower((*line)[0]) == "block" && (*line)[2] == "Q=") {
         const double block_scale = convert_to<double>((*line)[3]);
         if (flexiblesusy::is_equal(scale, block_scale))
            return true;
      }
   }

   return false;
}

std::string GM2_slha_io::to_lower(const std::string& str)
{
   std::string lower(str.size(), ' ');
   std::transform(str.begin(), str.end(), lower.begin(), ::tolower);
   return lower;
}

/**
 * Applies processor to each (key, value) pair of a SLHA block.
 * Non-data lines are ignored.
 *
 * @param block_name block name
 * @param processor tuple processor to be applied
 * @param scale (or 0 if scale should be ignored)
 */
void GM2_slha_io::read_block(const std::string& block_name,
                             const Tuple_processor& processor,
                             double scale) const
{
   SLHAea::Coll::const_iterator block =
      data.find(data.cbegin(), data.cend(), block_name);

   while (block != data.cend()) {
      if (at_scale(*block, scale)) {
         for (SLHAea::Block::const_iterator line = block->cbegin(),
                 end = block->cend(); line != end; ++line) {
            if (!line->is_data_line())
               continue;

            if (line->size() >= 2) {
               const int key = convert_to<int>((*line)[0]);
               const double value = convert_to<double>((*line)[1]);
               processor(key, value);
            }
         }
      }

      ++block;
      block = data.find(block, data.cend(), block_name);
   }
}

void GM2_slha_io::set_block(const std::ostringstream& lines, Position position)
{
   SLHAea::Block block;
   block.str(lines.str());
   data.erase(block.name());
   if (position == front)
      data.push_front(block);
   else
      data.push_back(block);
}

void GM2_slha_io::write_to_file(const std::string& file_name)
{
   std::ofstream ofs(file_name);
   write_to_stream(ofs);
}

void GM2_slha_io::write_to_stream(std::ostream& ostr)
{
   if (ostr.good())
      ostr << data;
   else
      ERROR("cannot write SLHA file");
}

/**
 * Fills a block entry with a value.  If the block or the entry do not
 * exist, the block / entry is created.
 *
 * @param block_name block name
 * @param entry number of the entry
 * @param value value
 * @param description comment
 */
void GM2_slha_io::fill_block_entry(const std::string& block_name,
                                   unsigned entry, double value,
                                   const std::string& description)
{
   std::ostringstream sstr;
   sstr << FORMAT_ELEMENT(entry, value, description);

   SLHAea::Coll::const_iterator block =
      data.find(data.cbegin(), data.cend(), block_name);

   if (block == data.cend()) {
      // create new block
      std::ostringstream block;
      block << "Block " << block_name << '\n'
            << sstr.str();
      set_block(block, GM2_slha_io::back);
   } else {
      data[block_name][SLHAea::to<double>(entry)] = sstr.str();
   }
}

/**
 * Fills a block entry with a string.  If the block or the entry do
 * not exist, the block / entry is created.
 *
 * @param block_name block name
 * @param entry number of the entry
 * @param description comment
 */
void GM2_slha_io::fill_block_entry(const std::string& block_name,
                                   unsigned entry,
                                   const std::string& description)
{
   std::ostringstream sstr;
   sstr << FORMAT_SPINFO(entry, description);

   SLHAea::Coll::const_iterator block =
      data.find(data.cbegin(), data.cend(), block_name);

   if (block == data.cend()) {
      // create new block
      std::ostringstream block;
      block << "Block " << block_name << '\n'
            << sstr.str();
      set_block(block, GM2_slha_io::back);
   } else {
      data[block_name][SLHAea::to<double>(entry)] = sstr.str();
   }
}

void fill_alpha_s(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   const double alpha_S = slha_io.read_entry("SMINPUTS", 3);

   if (!is_zero(alpha_S))
      model.set_g3(std::sqrt(4*M_PI*alpha_S));
}

void fill_soft_parameters_from_msoft(const GM2_slha_io& slha_io,
                                     MSSMNoFV_onshell& model, double scale)
{
   using namespace std::placeholders;

   GM2_slha_io::Tuple_processor processor
      = std::bind(process_msoft_tuple, std::ref(model), _1, _2);

   slha_io.read_block("MSOFT", processor, scale);
}

void fill_drbar_parameters(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   const double scale = slha_io.read_scale("HMIX");

   if (flexiblesusy::is_zero(scale)) {
      throw EInvalidInput("Could not determine renormalization scale"
                          " from HMIX block");
   }

   {
      Eigen::Matrix<double,3,3> Ae(Eigen::Matrix<double,3,3>::Zero());
      slha_io.read_block("AE", Ae, scale);
      model.set_Ae(Ae);
   }
   {
      Eigen::Matrix<double,3,3> Au(Eigen::Matrix<double,3,3>::Zero());
      slha_io.read_block("AU", Au, scale);
      model.set_Au(Au);
   }
   {
      Eigen::Matrix<double,3,3> Ad(Eigen::Matrix<double,3,3>::Zero());
      slha_io.read_block("AD", Ad, scale);
      model.set_Ad(Ad);
   }

   fill_soft_parameters_from_msoft(slha_io, model, scale);

   model.set_Mu(slha_io.read_entry("HMIX", 1, scale));

   const double tanb = slha_io.read_entry("HMIX", 2, scale);
   const double MA2_drbar = slha_io.read_entry("HMIX", 4, scale);
   const double sinb = tanb / std::sqrt(1 + tanb*tanb);
   const double cosb = 1.   / std::sqrt(1 + tanb*tanb);

   model.set_TB(tanb);
   model.set_BMu(MA2_drbar * sinb * cosb);
   model.set_scale(scale);
}

void fill_pole_masses_from_sminputs(
   const GM2_slha_io& slha_io, MSSMNoFV_onshell_physical& physical)
{
   using namespace std::placeholders;

   GM2_slha_io::Tuple_processor processor
      = std::bind(process_fermion_sminputs_tuple, std::ref(physical), _1, _2);

   slha_io.read_block("SMINPUTS", processor);
}

void fill_susy_masses_from_mass(
   const GM2_slha_io& slha_io, MSSMNoFV_onshell_physical& physical)
{
   using namespace std::placeholders;

   GM2_slha_io::Tuple_processor processor
      = std::bind(process_mass_tuple, std::ref(physical), _1, _2);

   slha_io.read_block("MASS", processor);
}

void fill_physical(const GM2_slha_io& slha_io, MSSMNoFV_onshell_physical& physical)
{
   // read all pole masses (includin MW) from SMINPUTS
   fill_pole_masses_from_sminputs(slha_io, physical);

   // if MW if given in MASS[24], prefer this value
   const double MW = slha_io.read_entry("MASS", 24);
   if (!is_zero(MW))
      physical.MVWm = MW;

   fill_susy_masses_from_mass(slha_io, physical);
}

void fill_pole_masses_from_sminputs_and_mass(
   const GM2_slha_io& slha_io, MSSMNoFV_onshell_physical& physical)
{
   MSSMNoFV_onshell_physical physical_hk(physical);
   fill_physical(slha_io, physical_hk);
   physical_hk.convert_to_hk();
   physical = physical_hk;
}

void fill_gm2_specific_alphas(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   const double alpha_MZ = std::abs(slha_io.read_entry("GM2CalcInput", 1));
   const double alpha_thompson = std::abs(slha_io.read_entry("GM2CalcInput", 2));

   if (alpha_MZ > std::numeric_limits<double>::epsilon())
      model.set_alpha_MZ(alpha_MZ);

   if (alpha_thompson > std::numeric_limits<double>::epsilon())
      model.set_alpha_thompson(alpha_thompson);
}

/**
 * Reads the GM2CalcInput block and fills the model parameter class.
 *
 * This function assumes that MW(pole) and MZ(pole) are non-zero.
 */
void fill_gm2_specific_onshell_parameters(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   using namespace std::placeholders;

   GM2_slha_io::Tuple_processor processor
      = std::bind(process_gm2calcinput_tuple, std::ref(model), _1, _2);

   slha_io.read_block("GM2CalcInput", processor);
}

/**
 * Reads model parameters in GM2Calc format from GM2CalcInput and
 * SMINPUTS blocks
 *
 * @param slha_io SLHA object
 * @param model model
 */
void fill_gm2calc(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   fill_pole_masses_from_sminputs(slha_io, model.get_physical());
   fill_alpha_s(slha_io, model);
   fill_gm2_specific_onshell_parameters(slha_io, model);
}

/**
 * Reads model parameters in SLHA format (from SLHA and GM2CalcInput
 * input blocks)
 *
 * @param slha_io SLHA object
 * @param model model
 */
void fill_slha(const GM2_slha_io& slha_io, MSSMNoFV_onshell& model)
{
   fill_pole_masses_from_sminputs_and_mass(slha_io, model.get_physical());
   fill_alpha_s(slha_io, model);
   fill_drbar_parameters(slha_io, model);
   fill_gm2_specific_alphas(slha_io, model);
}

/**
 * Reads configuration from GM2CalcConfig block
 *
 * @param slha_io SLHA object
 * @param config_options configuration settings
 */
void fill(const GM2_slha_io& slha_io, Config_options& config_options)
{
   using namespace std::placeholders;

   GM2_slha_io::Tuple_processor processor
      = std::bind(process_gm2calcconfig_tuple, std::ref(config_options), _1, _2);

   slha_io.read_block("GM2CalcConfig", processor);
}

namespace {
void process_gm2calcconfig_tuple(Config_options& config_options,
                                 int key, double value)
{
   switch (key) {
   case 0:
      config_options.output_format =
         static_cast<Config_options::E_output_format>(value);
      break;
   case 1:
      config_options.loop_order = value;
      break;
   case 2:
      config_options.tanb_resummation = value;
      break;
   case 3:
      config_options.force_output = value;
      break;
   case 4:
      config_options.verbose_output = value;
      break;
   default:
      WARNING("Unrecognized entry in block GM2CalcConfig: " << key);
      break;
   }
}

void process_gm2calcinput_tuple(MSSMNoFV_onshell& model,
                                             int key, double value)
{
   switch (key) {
   case 0: model.set_scale(value);          break;
   case 1: model.set_alpha_MZ(value);       break;
   case 2: model.set_alpha_thompson(value); break;
   case 3: {
      const double tanb = value;
      const double MW = model.get_MW();
      const double MZ = model.get_MZ();
      const double cW = MW/MZ;
      const double sW = std::sqrt(1. - cW*cW);
      const double vev = 2. * MW * sW / model.get_EL();
      const double sinb = tanb / std::sqrt(1 + tanb*tanb);
      const double cosb = 1.   / std::sqrt(1 + tanb*tanb);
      model.set_vd(vev * cosb);
      model.set_vu(vev * sinb);
      }
      break;
   case  4: model.set_Mu(    value); break;
   case  5: model.set_MassB( value); break;
   case  6: model.set_MassWB(value); break;
   case  7: model.set_MassG( value); break;
   case  8: model.set_MA0(   value); break;
   case  9: model.set_ml2(0, 0, signed_sqr(value)); break;
   case 10: model.set_ml2(1, 1, signed_sqr(value)); break;
   case 11: model.set_ml2(2, 2, signed_sqr(value)); break;
   case 12: model.set_me2(0, 0, signed_sqr(value)); break;
   case 13: model.set_me2(1, 1, signed_sqr(value)); break;
   case 14: model.set_me2(2, 2, signed_sqr(value)); break;
   case 15: model.set_mq2(0, 0, signed_sqr(value)); break;
   case 16: model.set_mq2(1, 1, signed_sqr(value)); break;
   case 17: model.set_mq2(2, 2, signed_sqr(value)); break;
   case 18: model.set_mu2(0, 0, signed_sqr(value)); break;
   case 19: model.set_mu2(1, 1, signed_sqr(value)); break;
   case 20: model.set_mu2(2, 2, signed_sqr(value)); break;
   case 21: model.set_md2(0, 0, signed_sqr(value)); break;
   case 22: model.set_md2(1, 1, signed_sqr(value)); break;
   case 23: model.set_md2(2, 2, signed_sqr(value)); break;
   case 24: model.set_Ae( 0, 0, value); break;
   case 25: model.set_Ae( 1, 1, value); break;
   case 26: model.set_Ae( 2, 2, value); break;
   case 27: model.set_Ad( 0, 0, value); break;
   case 28: model.set_Ad( 1, 1, value); break;
   case 29: model.set_Ad( 2, 2, value); break;
   case 30: model.set_Au( 0, 0, value); break;
   case 31: model.set_Au( 1, 1, value); break;
   case 32: model.set_Au( 2, 2, value); break;
   default:
      WARNING("Unrecognized entry in block GM2CalcInput: " << key);
      break;
   }
}

void process_fermion_sminputs_tuple(
   MSSMNoFV_onshell_physical& physical, int key, double value)
{
   switch (key) {
   case  1: /* alpha_em(MZ) */      break;
   case  2: /* G_F */               break;
   case  3: /* alpha_s(MZ) */       break;
   case  4: physical.MVZ = value;   break;
   case  5: physical.MFb = value;   break;
   case  6: physical.MFt = value;   break;
   case  7: physical.MFtau = value; break;
   case  8: /* nu_3 */              break;
   case  9: physical.MVWm = value;  break;
   case 11: physical.MFe = value;   break;
   case 12: /* nu_1 */              break;
   case 13: physical.MFm = value;   break;
   case 14: /* nu_2 */              break;
   case 21: physical.MFd = value;   break;
   case 23: physical.MFs = value;   break;
   case 22: physical.MFu = value;   break;
   case 24: physical.MFc = value;   break;
   default:
      WARNING("Unrecognized entry in block SMINPUTS: " << key);
      break;
   }
}

void process_msoft_tuple(
   MSSMNoFV_onshell& model, int key, double value)
{
   switch (key) {
   case 21: /* mHd2 */                            ; break;
   case 22: /* mHu2 */                            ; break;
   case 31: model.set_ml2(0, 0, signed_sqr(value)); break;
   case 32: model.set_ml2(1, 1, signed_sqr(value)); break;
   case 33: model.set_ml2(2, 2, signed_sqr(value)); break;
   case 34: model.set_me2(0, 0, signed_sqr(value)); break;
   case 35: model.set_me2(1, 1, signed_sqr(value)); break;
   case 36: model.set_me2(2, 2, signed_sqr(value)); break;
   case 41: model.set_mq2(0, 0, signed_sqr(value)); break;
   case 42: model.set_mq2(1, 1, signed_sqr(value)); break;
   case 43: model.set_mq2(2, 2, signed_sqr(value)); break;
   case 44: model.set_mu2(0, 0, signed_sqr(value)); break;
   case 45: model.set_mu2(1, 1, signed_sqr(value)); break;
   case 46: model.set_mu2(2, 2, signed_sqr(value)); break;
   case 47: model.set_md2(0, 0, signed_sqr(value)); break;
   case 48: model.set_md2(1, 1, signed_sqr(value)); break;
   case 49: model.set_md2(2, 2, signed_sqr(value)); break;
   case  1: model.set_MassB(               value) ; break;
   case  2: model.set_MassWB(              value) ; break;
   case  3: model.set_MassG(               value) ; break;
   default:
      WARNING("Unrecognized entry in block MSOFT: " << key);
      break;
   }
}

void process_mass_tuple(
   MSSMNoFV_onshell_physical& physical, int key, double value)
{
   switch (key) {
   case 1000012: physical.MSveL = value;    break;
   case 1000014: physical.MSvmL = value;    break;
   case 1000016: physical.MSvtL = value;    break;
   case 1000001: physical.MSd(0) = value;   break;
   case 2000001: physical.MSd(1) = value;   break;
   case 1000002: physical.MSu(0) = value;   break;
   case 2000002: physical.MSu(1) = value;   break;
   case 1000011: physical.MSe(0) = value;   break;
   case 2000011: physical.MSe(1) = value;   break;
   case 1000013: physical.MSm(0) = value;   break;
   case 2000013: physical.MSm(1) = value;   break;
   case 1000015: physical.MStau(0) = value; break;
   case 2000015: physical.MStau(1) = value; break;
   case 1000003: physical.MSs(0) = value;   break;
   case 2000003: physical.MSs(1) = value;   break;
   case 1000004: physical.MSc(0) = value;   break;
   case 2000004: physical.MSc(1) = value;   break;
   case 1000005: physical.MSb(0) = value;   break;
   case 2000005: physical.MSb(1) = value;   break;
   case 1000006: physical.MSt(0) = value;   break;
   case 2000006: physical.MSt(1) = value;   break;
   case 24     : /* MW */                   break;
   case 25     : physical.Mhh(0) = value;   break;
   case 35     : physical.Mhh(1) = value;   break;
   case 36     : physical.MAh(1) = value;   break;
   case 37     : physical.MHpm(1) = value;  break;
   case 1000021: /* gluino */               break;
   case 1000022: physical.MChi(0) = value;  break;
   case 1000023: physical.MChi(1) = value;  break;
   case 1000025: physical.MChi(2) = value;  break;
   case 1000035: physical.MChi(3) = value;  break;
   case 1000024: physical.MCha(0) = value;  break;
   case 1000037: physical.MCha(1) = value;  break;
   default:
      break;
   }
}
} // anonymous namespace

} // namespace gm2calc
