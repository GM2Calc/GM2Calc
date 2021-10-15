#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"

#include "gm2_config_options.hpp"
#include "gm2_numerics.hpp"
#include "gm2_slha_io.hpp"

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

   CHECK_THROWS(slha.fill(config));

   CHECK(config.output_format         == default_config.output_format);
   CHECK(config.loop_order            == default_config.loop_order);
   CHECK(config.tanb_resummation      == default_config.tanb_resummation);
   CHECK(config.force_output          == default_config.force_output);
   CHECK(config.verbose_output        == default_config.verbose_output);
   CHECK(config.calculate_uncertainty == default_config.calculate_uncertainty);
}


TEST_CASE("read_matrix_dense")
{
   const double eps = std::numeric_limits<double>::epsilon();

   char const * const slha_input = R"(
BLOCK MAT Q= 100
    1   1  1
    1   2  2
    2   1  3
    2   2  4
    3   1  5
    3   2  6
)";

   Eigen::Matrix<double,3,2> m;
   m.setZero();

   std::istringstream stream(slha_input);
   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.read_block("MAT", m);

   CHECK_CLOSE(m(0,0), 1.0, eps);
   CHECK_CLOSE(m(0,1), 2.0, eps);
   CHECK_CLOSE(m(1,0), 3.0, eps);
   CHECK_CLOSE(m(1,1), 4.0, eps);
   CHECK_CLOSE(m(2,0), 5.0, eps);
   CHECK_CLOSE(m(2,1), 6.0, eps);
}


TEST_CASE("read_vector_dense")
{
   const double eps = std::numeric_limits<double>::epsilon();

   char const * const slha_input = R"(
BLOCK VEC Q= 100
    1   1
    2   2
    3   3
)";

   Eigen::Matrix<double,3,1> v;
   v.setZero();

   std::istringstream stream(slha_input);
   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.read_block("VEC", v);

   CHECK_CLOSE(v(0), 1.0, eps);
   CHECK_CLOSE(v(1), 2.0, eps);
   CHECK_CLOSE(v(2), 3.0, eps);
}


TEST_CASE("fill_block_entry")
{
   gm2calc::GM2_slha_io slha;

   slha.fill_block_entry("MAT", 1, 2.0, "description");
   slha.fill_block_entry("INFO", 1, "description");

   slha.write_to_stream(std::cout);

   slha.fill_block_entry("MAT", 2, 3.0, "description");
   slha.fill_block_entry("INFO", 2, "description");

   slha.write_to_stream(std::cout);
}


TEST_CASE("read_scale")
{
   const double eps = std::numeric_limits<double>::epsilon();

   char const * const slha_input = R"(
BLOCK VEC Q= 50
    1   1
    2   2
    3   3
BLOCK VEC Q= 100
    1   1
    2   2
    3   3
)";

   std::istringstream stream(slha_input);
   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   const double scale = slha.read_scale("VEC");
   CHECK_CLOSE(scale, 100.0, eps);
}


TEST_CASE("read_scale_no_scale")
{
   const double eps = std::numeric_limits<double>::epsilon();

   char const * const slha_input = R"(
BLOCK VEC
    1   1
    2   2
    3   3
)";

   std::istringstream stream(slha_input);
   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   const double scale = slha.read_scale("VEC");
   CHECK_CLOSE(scale, 0.0, eps);
}
