#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_slha_io.hpp"
#include "gm2_config_options.hpp"

#include <sstream>

#define CHECK_CLOSE(a,b,eps)                            \
   do {                                                 \
      CHECK((a) == doctest::Approx(b).epsilon(eps));    \
   } while (0)


TEST_CASE("fill_valid_config")
{
   gm2calc::Config_options config;

   char const * const slha_input = R"(
Block GM2CalcConfig
     0     1     # output format (0, 1, 2, 3, 4)
     1     1     # loop order (0, 1 or 2)
     2     0     # disable/enable tan(beta) resummation (0 or 1)
     3     1     # force output (0 or 1)
     4     1     # verbose output (0 or 1)
     5     1     # calculate uncertainty
)";

   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill(config);

   CHECK(config.output_format == 1);
   CHECK(config.loop_order == 1);
   CHECK(config.tanb_resummation == false);
   CHECK(config.force_output == true);
   CHECK(config.verbose_output == true);
   CHECK(config.calculate_uncertainty == true);
}

// checks that default configuration values are used if invalid
// configuration values are passed
TEST_CASE("fill_invalid_config")
{
   const gm2calc::Config_options default_config;
   gm2calc::Config_options config;

   char const * const slha_input = R"(
Block GM2CalcConfig
     0     9     # output format (0, 1, 2, 3, 4)
     1     9     # loop order (0, 1 or 2)
     2     2     # disable/enable tan(beta) resummation (0 or 1)
     3     2     # force output (0 or 1)
     4     2     # verbose output (0 or 1)
     5     2     # calculate uncertainty (0 or 1)
)";

   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill(config);

   CHECK(config.output_format         == default_config.output_format);
   CHECK(config.loop_order            == default_config.loop_order);
   CHECK(config.tanb_resummation      == default_config.tanb_resummation);
   CHECK(config.force_output          == default_config.force_output);
   CHECK(config.verbose_output        == default_config.verbose_output);
   CHECK(config.calculate_uncertainty == default_config.calculate_uncertainty);
}
