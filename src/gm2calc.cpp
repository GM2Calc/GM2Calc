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

#include "gm2_1loop.hpp"
#include "gm2_2loop.hpp"
#include "gm2_error.hpp"
#include "gm2_uncertainty.hpp"
#include "config.h"

#include "gm2_slha_io.hpp"
#include "MSSMNoFV_onshell.hpp"

#include <iostream>
#include <limits>

#define ERROR(message) std::cerr << "Error: " << message << '\n';
#define WARNING(message) std::cerr << "Warning: " << message << '\n';

using namespace flexiblesusy;

/**
 * @class Gm2_cmd_line_options
 * @brief command line options for GM2Calc
 */
struct Gm2_cmd_line_options {
   enum E_input_type { SLHA, GM2Calc };

   std::string input_source; ///< input source (file name or `-' for stdin)
   E_input_type input_type;  ///< input format (SLHA or GM2Calc)

   static bool starts_with(const std::string& str, const std::string& prefix) {
      return !str.compare(0, prefix.size(), prefix);
   }
};

void print_usage(const char* program_name)
{
   std::cout <<
      "Usage: " << program_name << " [options]\n"
      "Options:\n"
      "  --slha-input-file=<source>      SLHA input source (file name or - for stdin)\n"
      "  --gm2calc-input-file=<source>   GM2Calc input source (file name or - for stdin)\n"
      "  --help,-h                       print this help message\n"
      "  --version,-v                    print version number"
             << std::endl;
}

Gm2_cmd_line_options get_cmd_line_options(int argc, const char* argv[])
{
   Gm2_cmd_line_options options;

   for (int i = 1; i < argc; ++i) {
      const std::string option(argv[i]);

      if (Gm2_cmd_line_options::starts_with(option, "--slha-input-file=")) {
         options.input_source = option.substr(18);
         options.input_type = Gm2_cmd_line_options::SLHA;
         continue;
      }

      if (Gm2_cmd_line_options::starts_with(option, "--gm2calc-input-file=")) {
         options.input_source = option.substr(21);
         options.input_type = Gm2_cmd_line_options::GM2Calc;
         continue;
      }

      if (option == "--help" || option == "-h") {
         print_usage(argv[0]);
         exit(EXIT_SUCCESS);
      }

      if (option == "--version" || option == "-v") {
         std::cout << GM2CALC_VERSION << '\n';
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }

   return options;
}

/**
 * Set the config options to default values, depending on the input
 * parameter set (chosen by the user).
 *
 * If SLHA input format has been selected, the default output format
 * will be SLHA format.  By default \f$a_\mu\f$ will be written to the
 * SPheno block SPhenoLowEnergy[21].
 *
 * If GM2Calc input format has been chosen, the default values set in
 * \a Config_options are used.
 */
void set_to_default(gm2calc::Config_options& config_options,
                    const Gm2_cmd_line_options& options)
{
   switch (options.input_type) {
   case Gm2_cmd_line_options::SLHA:
      config_options.output_format = gm2calc::Config_options::SPheno;
      break;
   case Gm2_cmd_line_options::GM2Calc:
      break;
   default:
      throw gm2calc::ESetupError("Unknown input option");
      break;
   }
}

/**
 * Setup the model parameters consistently, depending on the input
 * parameter set (chosen by the user).
 *
 * If the user has chosen an SLHA-compliant input, the input
 * parameters (given in \a slha_io) are converted to the g-2 specific
 * on-shell scheme.
 *
 * If the user has chosen the GM2Calc specific input, the input
 * parameters are already given in the on-shell scheme.  In this case
 * the SUSY mass spectrum is calculated from the given input
 * parameters.
 *
 * @param model model parameters (and particle masses)
 * @param slha_io object with numerical values of input parameters
 * @param options command line options (defines meaning of input parameter set)
 * @param config_options configuration options for the calculation
 */
void setup_model(gm2calc::MSSMNoFV_onshell& model,
                 const gm2calc::GM2_slha_io& slha_io,
                 const Gm2_cmd_line_options& options,
                 const gm2calc::Config_options& config_options)
{
   model.set_verbose_output(config_options.verbose_output);

   switch (options.input_type) {
   case Gm2_cmd_line_options::SLHA:
      // determine on-shell model parameters from an SLHA parameter
      // set
      fill_slha(slha_io, model);
      model.convert_to_onshell();
      break;
   case Gm2_cmd_line_options::GM2Calc:
      // on-shell parameters are directly given, calculate mass
      // spectrum
      fill_gm2calc(slha_io, model);
      model.calculate_masses();
      break;
   default:
      throw gm2calc::ESetupError("Unknown input option");
      break;
   }

   if (model.do_verbose_output())
      std::cout << model << '\n';
}

/**
 * Prints detailed a_mu calculation (1-loop w/ and w/o tan(beta)
 * resummation, 2-loop, and different contributions).
 *
 * @param model model object (contains parameters)
 */
void print_amu_detailed(
   const gm2calc::MSSMNoFV_onshell& model)
{
#define FORMAT_AMU(amu) boost::format("% 16.14e") % (amu)
#define FORMAT_DEL(amu) boost::format("%16.14e") % (amu)
#define FORMAT_PCT(pct) boost::format("%2.1f") % (pct)

   std::string error_str;
   if (model.get_problems().have_problem()) {
      error_str = model.get_problems().get_problems()
         + " (with tan(beta) resummation)\n\n";
   }

   const double amu_1l = gm2calc::calculate_amu_1loop(model);
   const double amu_2l_photonic_chipm = gm2calc::amuChipmPhotonic(model);
   const double amu_2l_photonic_chi0 = gm2calc::amuChi0Photonic(model);
   const double amu_2l_a_sfermion = gm2calc::amu2LaSferm(model);
   const double amu_2l_a_cha = gm2calc::amu2LaCha(model);
   const double amu_2l_ferm_sferm_approx = gm2calc::amu2LFSfapprox(model);
   const double amu_2l = gm2calc::calculate_amu_2loop(model);
   const double amu_2l_uncertainty = gm2calc::calculate_uncertainty_amu_2loop(model);
   const double tan_beta_cor = gm2calc::tan_beta_cor(model);

   // no tan(beta) resummation
   double amu_1l_non_tan_beta_resummed = 0.;
   double amu_2l_non_tan_beta_resummed = 0.;
   std::string error_str_non_tan_beta_resummation;

   try {
      // w/o tan(beta) resummation, allow throwing exceptions
      gm2calc::MSSMNoFV_onshell model_except(model);
      model_except.do_force_output(false);
      amu_1l_non_tan_beta_resummed =
         gm2calc::calculate_amu_1loop_non_tan_beta_resummed(model_except);
      amu_2l_non_tan_beta_resummed =
         gm2calc::calculate_amu_2loop_non_tan_beta_resummed(model_except);
   } catch (const gm2calc::Error& error) {
      error_str_non_tan_beta_resummation = " (" + error.what() + ")";
      // try to redo calculation w/o throwing an exception
      gm2calc::MSSMNoFV_onshell model_no_except(model);
      model_no_except.do_force_output(true);
      amu_1l_non_tan_beta_resummed =
         gm2calc::calculate_amu_1loop_non_tan_beta_resummed(model_no_except);
      amu_2l_non_tan_beta_resummed =
         gm2calc::calculate_amu_2loop_non_tan_beta_resummed(model_no_except);
   }

   const double amu_2l_tanb_approx =
      + (tan_beta_cor - 1.) * amu_1l_non_tan_beta_resummed;

   const double amu_best = amu_1l + amu_2l;

   std::cout <<
      "===============================================================================\n"
      "   amu (1-loop + 2-loop best) = " << FORMAT_AMU(amu_best) << ' ' <<
      "+- " << FORMAT_DEL(amu_2l_uncertainty) << '\n' <<
      "===============================================================================\n"
      "\n" <<
      error_str <<
      "==============================\n"
      "   amu (1-loop) corrections\n"
      "==============================\n"
      "\n"
      "full 1L with tan(beta) resummation:\n"
      "   chi^0     " << FORMAT_AMU(gm2calc::amuChi0(model)) << '\n' <<
      "   chi^+-    " << FORMAT_AMU(gm2calc::amuChipm(model)) << '\n' <<
      "   -------------------------------\n"
      "   sum       " << FORMAT_AMU(amu_1l) <<
      " (" << FORMAT_PCT(100. * amu_1l / amu_best) << "% of full 1L + 2L result)\n"
      "\n"
      "full 1L without tan(beta) resummation:\n"
      "             " << FORMAT_AMU(amu_1l_non_tan_beta_resummed) <<
      error_str_non_tan_beta_resummation << '\n' <<
      "\n"
      "1L approximation with tan(beta) resummation:\n"
      "   W-H-nu    " << FORMAT_AMU(gm2calc::amuWHnu(model) * tan_beta_cor) << '\n' <<
      "   W-H-muL   " << FORMAT_AMU(gm2calc::amuWHmuL(model) * tan_beta_cor) << '\n' <<
      "   B-H-muL   " << FORMAT_AMU(gm2calc::amuBHmuL(model) * tan_beta_cor) << '\n' <<
      "   B-H-muR   " << FORMAT_AMU(gm2calc::amuBHmuR(model) * tan_beta_cor) << '\n' <<
      "   B-muL-muR " << FORMAT_AMU(gm2calc::amuBmuLmuR(model) * tan_beta_cor) << '\n' <<
      "   -------------------------------\n"
      "   sum       " << FORMAT_AMU(gm2calc::amu1Lapprox(model)) << '\n' <<
      "\n"
      "==============================\n"
      "   amu (2-loop) corrections\n"
      "==============================\n"
      "\n"
      "2L best with tan(beta) resummation:\n"
      "             " << FORMAT_AMU(amu_2l) <<
      " (" << FORMAT_PCT(100. * amu_2l / amu_best) << "% of full 1L + 2L result)\n"
      "\n"
      "2L best without tan(beta) resummation:\n"
      "             " << FORMAT_AMU(amu_2l_non_tan_beta_resummed) <<
      error_str_non_tan_beta_resummation << '\n' <<
      "\n"
      "photonic with tan(beta) resummation:\n"
      "   chi^0     " << FORMAT_AMU(amu_2l_photonic_chi0) << '\n' <<
      "   chi^+-    " << FORMAT_AMU(amu_2l_photonic_chipm) << '\n' <<
      "   -------------------------------\n"
      "   sum       " << FORMAT_AMU(amu_2l_photonic_chipm
                                    + amu_2l_photonic_chi0) <<
      " (" << FORMAT_PCT(100. * (amu_2l_photonic_chipm + amu_2l_photonic_chi0)
                         / amu_best) << "% of full 1L + 2L result)\n"
      "\n"
      "fermion/sfermion approximation with tan(beta) resummation:\n"
      "   W-H-nu    " << FORMAT_AMU(gm2calc::amuWHnu2L(model) * tan_beta_cor) << '\n' <<
      "   W-H-muL   " << FORMAT_AMU(gm2calc::amuWHmuL2L(model) * tan_beta_cor) << '\n' <<
      "   B-H-muL   " << FORMAT_AMU(gm2calc::amuBHmuL2L(model) * tan_beta_cor) << '\n' <<
      "   B-H-muR   " << FORMAT_AMU(gm2calc::amuBHmuR2L(model)* tan_beta_cor) << '\n' <<
      "   B-muL-muR " << FORMAT_AMU(gm2calc::amuBmuLmuR2L(model) * tan_beta_cor) << '\n' <<
      "   -------------------------------\n"
      "   sum       " << FORMAT_AMU(amu_2l_ferm_sferm_approx) <<
      " (" << FORMAT_PCT(100. * amu_2l_ferm_sferm_approx / amu_best) << "% of full 1L + 2L result)\n"
      "\n"
      "2L(a) (1L insertions into 1L SM diagram) with tan(beta) resummation:\n"
      "   sfermion  " << FORMAT_AMU(amu_2l_a_sfermion) << '\n' <<
      "   cha^+-    " << FORMAT_AMU(amu_2l_a_cha) << '\n' <<
      "   -------------------------------\n"
      "   sum       " << FORMAT_AMU(amu_2l_a_sfermion + amu_2l_a_cha) <<
      " (" << FORMAT_PCT(100. * (amu_2l_a_sfermion + amu_2l_a_cha) / amu_best) << "% of full 1L + 2L result)\n"
      "\n"
      "tan(beta) correction:\n"
      "   amu(1L) * (1 / (1 + Delta_mu) - 1) = " << FORMAT_AMU(amu_2l_tanb_approx) <<
      " (" << FORMAT_PCT(100. * amu_2l_tanb_approx / amu_1l_non_tan_beta_resummed) << "%)\n"
      ;

#undef FORMAT_AMU
#undef FORMAT_PCT
}

/**
 * Calculates a_mu for a given set of configuration options (loop
 * order, tan(beta) resummation.
 *
 * @param model model parameters
 * @param config_options configuration options
 *
 * @return a_mu
 */
double calculate_amu(const gm2calc::MSSMNoFV_onshell& model,
                     const gm2calc::Config_options& config_options)
{
   double result = std::numeric_limits<double>::signaling_NaN();

   switch (config_options.loop_order) {
   case 0:
      result = 0.;
      break;
   case 1:
      if (config_options.tanb_resummation)
         result = gm2calc::calculate_amu_1loop(model);
      else
         result = gm2calc::calculate_amu_1loop_non_tan_beta_resummed(model);
      break;
   case 2:
      if (config_options.tanb_resummation)
         result = gm2calc::calculate_amu_1loop(model)
            + gm2calc::calculate_amu_2loop(model);
      else
         result = gm2calc::calculate_amu_1loop_non_tan_beta_resummed(model)
            + gm2calc::calculate_amu_2loop_non_tan_beta_resummed(model);
      break;
   default:
      ERROR("loop orders > 2 not supported!");
      break;
   }

   return result;
}

/**
 * Calculates amu and prints it to std::cout.  The output format
 * depends on the \a config_options .
 *
 * @param model model
 * @param slha_io SLHA object
 * @param config_options configuration options
 */
void print_amu(const gm2calc::MSSMNoFV_onshell& model,
               gm2calc::GM2_slha_io& slha_io,
               const gm2calc::Config_options& config_options)
{
   switch (config_options.output_format) {
   case gm2calc::Config_options::Minimal:
      std::cout << std::setprecision(std::numeric_limits<double>::digits10)
                << std::scientific
                << (!config_options.calculate_uncertainty ?
                    calculate_amu(model, config_options) :
                    calculate_uncertainty_amu_2loop(model))
                << '\n';
      break;
   case gm2calc::Config_options::Detailed:
      print_amu_detailed(model);
      break;
   case gm2calc::Config_options::NMSSMTools:
      slha_io.fill_block_entry("LOWEN", 6,
                               calculate_amu(model, config_options),
                               "Delta(g-2)_muon/2");
      if (config_options.calculate_uncertainty) {
         slha_io.fill_block_entry("GM2CalcOutput", 1,
                                  calculate_uncertainty_amu_2loop(model),
                                  "uncertainty of a_mu(2-loop, tan(beta) resummation)");
      }
      slha_io.write_to_stream(std::cout);
      break;
   case gm2calc::Config_options::SPheno:
      slha_io.fill_block_entry("SPhenoLowEnergy", 21,
                               calculate_amu(model, config_options),
                               "Delta(g-2)_muon/2");
      if (config_options.calculate_uncertainty) {
         slha_io.fill_block_entry("GM2CalcOutput", 1,
                                  calculate_uncertainty_amu_2loop(model),
                                  "uncertainty of a_mu(2-loop, tan(beta) resummation)");
      }
      slha_io.write_to_stream(std::cout);
      break;
   default:
      ERROR("Unknown output format: " << config_options.output_format);
      break;
   }
}

/**
 * Prints output if an error has occured.
 *
 * @param error error object
 * @param slha_io SLHA object
 * @param config_options configuration options
 */
void print_error(const gm2calc::Error& error,
                 gm2calc::GM2_slha_io& slha_io,
                 const gm2calc::Config_options& config_options)
{
   ERROR(error.what());

   switch (config_options.output_format) {
   case gm2calc::Config_options::NMSSMTools:
   case gm2calc::Config_options::SPheno:
      // print SPINFO block with error description
      slha_io.fill_block_entry("SPINFO", 1, "GM2Calc");
      slha_io.fill_block_entry("SPINFO", 2, GM2CALC_VERSION);
      slha_io.fill_block_entry("SPINFO", 4, error.what());
      slha_io.write_to_stream(std::cout);
      break;
   default:
      break;
   }
}

/**
 * Adds SPINFO block with warning message to SLHA output if there is a
 * warning.
 *
 * @param model model
 * @param slha_io SLHA object
 * @param config_options configuration options
 */
void print_warnings(const gm2calc::MSSMNoFV_onshell& model,
                    gm2calc::GM2_slha_io& slha_io,
                    const gm2calc::Config_options& config_options)
{
   if (model.get_problems().have_problem() ||
       model.get_problems().have_warning())
      WARNING(model.get_problems());

   if (model.get_problems().have_warning()) {
      switch (config_options.output_format) {
      case gm2calc::Config_options::NMSSMTools:
      case gm2calc::Config_options::SPheno:
         // print SPINFO block with warning description
         slha_io.fill_block_entry("SPINFO", 1, "GM2Calc");
         slha_io.fill_block_entry("SPINFO", 2, GM2CALC_VERSION);
         slha_io.fill_block_entry("SPINFO", 3, model.get_problems().get_warnings());
         break;
      default:
         break;
      }
   }
}

int main(int argc, const char* argv[])
{
   Gm2_cmd_line_options options(get_cmd_line_options(argc, argv));

   if (options.input_source.empty()) {
      ERROR("No input source given!\n"
            "   Please provide an SLHA input via the option --slha-input-file=\n"
            "   or a GM2Calc input via the option --gm2calc-input-file=");
      return EXIT_FAILURE;
   }

   gm2calc::MSSMNoFV_onshell model;
   gm2calc::GM2_slha_io slha_io;
   gm2calc::Config_options config_options;
   int exit_code = EXIT_SUCCESS;

   try {
      set_to_default(config_options, options);
      slha_io.read_from_source(options.input_source);
      fill(slha_io, config_options);
      model.do_force_output(config_options.force_output);
      setup_model(model, slha_io, options, config_options);
      print_warnings(model, slha_io, config_options);
      print_amu(model, slha_io, config_options);
   } catch (const gm2calc::Error& error) {
      print_error(error, slha_io, config_options);
      exit_code = EXIT_FAILURE;
   }

   if (model.get_problems().have_problem())
      exit_code = EXIT_FAILURE;

   return exit_code;
}
