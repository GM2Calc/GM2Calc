#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2calc/THDM.hpp"
#include "gm2_slha_io.hpp"
#include <sstream>

TEST_CASE("fill_gauge_basis")
{
   char const * const slha_input = R"(
Block MINPAR                # model parameters
     3     3                # tan(beta)
    11     0.7              # lambda_1
    12     0.6              # lambda_2
    13     0.5              # lambda_3
    14     0.4              # lambda_4
    15     0.3              # lambda_5
    16     0.2              # lambda_6
    17     0.1              # lambda_7
    18     40000            # m_{12}^2
    21     0.01             # zeta_u
    22     0.02             # zeta_d
    23     0.03             # zeta_l
    24     2                # Yukawa type (0, 1, 2, 3, 4, 5)
Block GM2CalcTHDMXuInput
    1 1   1                    # Re(Xu(1,1))
    1 2   2                    # Re(Xu(1,2))
    1 3   3                    # Re(Xu(1,3))
    2 1   4                    # Re(Xu(2,1))
    2 2   5                    # Re(Xu(2,2))
    2 3   6                    # Re(Xu(2,3))
    3 1   7                    # Re(Xu(3,1))
    3 2   8                    # Re(Xu(3,2))
    3 3   9                    # Re(Xu(3,3))
Block GM2CalcTHDMXdInput
    1 1   11                   # Re(Xd(1,1))
    1 2   12                   # Re(Xd(1,2))
    1 3   13                   # Re(Xd(1,3))
    2 1   14                   # Re(Xd(2,1))
    2 2   15                   # Re(Xd(2,2))
    2 3   16                   # Re(Xd(2,3))
    3 1   17                   # Re(Xd(3,1))
    3 2   18                   # Re(Xd(3,2))
    3 3   19                   # Re(Xd(3,3))
Block GM2CalcTHDMXlInput
    1 1   21                   # Re(Xl(1,1))
    1 2   22                   # Re(Xl(1,2))
    1 3   23                   # Re(Xl(1,3))
    2 1   24                   # Re(Xl(2,1))
    2 2   25                   # Re(Xl(2,2))
    2 3   26                   # Re(Xl(2,3))
    3 1   27                   # Re(Xl(3,1))
    3 2   28                   # Re(Xl(3,2))
    3 3   29                   # Re(Xl(3,3))
Block GM2CalcTHDMImXuInput
    1 1   31                    # Im(Xu(1,1))
    1 2   32                    # Im(Xu(1,2))
    1 3   33                    # Im(Xu(1,3))
    2 1   34                    # Im(Xu(2,1))
    2 2   35                    # Im(Xu(2,2))
    2 3   36                    # Im(Xu(2,3))
    3 1   37                    # Im(Xu(3,1))
    3 2   38                    # Im(Xu(3,2))
    3 3   39                    # Im(Xu(3,3))
Block GM2CalcTHDMImXdInput
    1 1   41                    # Im(Xd(1,1))
    1 2   42                    # Im(Xd(1,2))
    1 3   43                    # Im(Xd(1,3))
    2 1   44                    # Im(Xd(2,1))
    2 2   45                    # Im(Xd(2,2))
    2 3   46                    # Im(Xd(2,3))
    3 1   47                    # Im(Xd(3,1))
    3 2   48                    # Im(Xd(3,2))
    3 3   49                    # Im(Xd(3,3))
Block GM2CalcTHDMImXlInput
    1 1   51                    # Im(Xl(1,1))
    1 2   52                    # Im(Xl(1,2))
    1 3   53                    # Im(Xl(1,3))
    2 1   54                    # Im(Xl(2,1))
    2 2   55                    # Im(Xl(2,2))
    2 3   56                    # Im(Xl(2,3))
    3 1   57                    # Im(Xl(3,1))
    3 2   58                    # Im(Xl(3,2))
    3 3   59                    # Im(Xl(3,3))
)";

   gm2calc::thdm::Gauge_basis basis;
   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill(basis);

   CHECK(basis.tan_beta == 3.0);
   CHECK(basis.lambda(0) == 0.7);
   CHECK(basis.lambda(1) == 0.6);
   CHECK(basis.lambda(2) == 0.5);
   CHECK(basis.lambda(3) == 0.4);
   CHECK(basis.lambda(4) == 0.3);
   CHECK(basis.lambda(5) == 0.2);
   CHECK(basis.lambda(6) == 0.1);
   CHECK(basis.m122 == 40000.0);
   CHECK(basis.zeta_u == 0.01);
   CHECK(basis.zeta_d == 0.02);
   CHECK(basis.zeta_l == 0.03);
   CHECK(basis.yukawa_type == gm2calc::thdm::Yukawa_type::type_2);

   for (int i = 0, count = 1; i < 3; i++) {
      for (int k = 0; k < 3; k++, count++) {
         CHECK(std::real(basis.Xu(i,k)) == count);
         CHECK(std::real(basis.Xd(i,k)) == (count + 10));
         CHECK(std::real(basis.Xl(i,k)) == (count + 20));
         CHECK(std::imag(basis.Xu(i,k)) == (count + 30));
         CHECK(std::imag(basis.Xd(i,k)) == (count + 40));
         CHECK(std::imag(basis.Xl(i,k)) == (count + 50));
      }
   }
}


TEST_CASE("fill_mass_basis")
{
   char const * const slha_input = R"(
Block MINPAR                # model parameters
     3     3                # tan(beta)
    16     0.2              # lambda_6
    17     0.1              # lambda_7
    18     40000            # m_{12}^2
    20     0.995            # sin(beta - alpha)
    21     0.01             # zeta_u
    22     0.02             # zeta_d
    23     0.03             # zeta_l
    24     2                # Yukawa type (0, 1, 2, 3, 4, 5)
Block MASS                  # Higgs masses
    25     100              # mh, lightest CP-even Higgs
    35     400              # mH, heaviest CP-even Higgs
    36     420              # mA, CP-odd Higgs
    37     440              # mH+, charged Higgs
Block GM2CalcTHDMXuInput
    1 1   1                    # Re(Xu(1,1))
    1 2   2                    # Re(Xu(1,2))
    1 3   3                    # Re(Xu(1,3))
    2 1   4                    # Re(Xu(2,1))
    2 2   5                    # Re(Xu(2,2))
    2 3   6                    # Re(Xu(2,3))
    3 1   7                    # Re(Xu(3,1))
    3 2   8                    # Re(Xu(3,2))
    3 3   9                    # Re(Xu(3,3))
Block GM2CalcTHDMXdInput
    1 1   11                   # Re(Xd(1,1))
    1 2   12                   # Re(Xd(1,2))
    1 3   13                   # Re(Xd(1,3))
    2 1   14                   # Re(Xd(2,1))
    2 2   15                   # Re(Xd(2,2))
    2 3   16                   # Re(Xd(2,3))
    3 1   17                   # Re(Xd(3,1))
    3 2   18                   # Re(Xd(3,2))
    3 3   19                   # Re(Xd(3,3))
Block GM2CalcTHDMXlInput
    1 1   21                   # Re(Xl(1,1))
    1 2   22                   # Re(Xl(1,2))
    1 3   23                   # Re(Xl(1,3))
    2 1   24                   # Re(Xl(2,1))
    2 2   25                   # Re(Xl(2,2))
    2 3   26                   # Re(Xl(2,3))
    3 1   27                   # Re(Xl(3,1))
    3 2   28                   # Re(Xl(3,2))
    3 3   29                   # Re(Xl(3,3))
Block GM2CalcTHDMImXuInput
    1 1   31                    # Im(Xu(1,1))
    1 2   32                    # Im(Xu(1,2))
    1 3   33                    # Im(Xu(1,3))
    2 1   34                    # Im(Xu(2,1))
    2 2   35                    # Im(Xu(2,2))
    2 3   36                    # Im(Xu(2,3))
    3 1   37                    # Im(Xu(3,1))
    3 2   38                    # Im(Xu(3,2))
    3 3   39                    # Im(Xu(3,3))
Block GM2CalcTHDMImXdInput
    1 1   41                    # Im(Xd(1,1))
    1 2   42                    # Im(Xd(1,2))
    1 3   43                    # Im(Xd(1,3))
    2 1   44                    # Im(Xd(2,1))
    2 2   45                    # Im(Xd(2,2))
    2 3   46                    # Im(Xd(2,3))
    3 1   47                    # Im(Xd(3,1))
    3 2   48                    # Im(Xd(3,2))
    3 3   49                    # Im(Xd(3,3))
Block GM2CalcTHDMImXlInput
    1 1   51                    # Im(Xl(1,1))
    1 2   52                    # Im(Xl(1,2))
    1 3   53                    # Im(Xl(1,3))
    2 1   54                    # Im(Xl(2,1))
    2 2   55                    # Im(Xl(2,2))
    2 3   56                    # Im(Xl(2,3))
    3 1   57                    # Im(Xl(3,1))
    3 2   58                    # Im(Xl(3,2))
    3 3   59                    # Im(Xl(3,3))
)";

   gm2calc::thdm::Mass_basis basis;
   std::istringstream stream(slha_input);

   gm2calc::GM2_slha_io slha;
   slha.read_from_stream(stream);
   slha.fill(basis);

   CHECK(basis.tan_beta == 3.0);
   CHECK(basis.mh == 100.0);
   CHECK(basis.mH == 400.0);
   CHECK(basis.mA == 420.0);
   CHECK(basis.mHp == 440.0);
   CHECK(basis.lambda6 == 0.2);
   CHECK(basis.lambda7 == 0.1);
   CHECK(basis.sin_beta_minus_alpha == 0.995);
   CHECK(basis.m122 == 40000.0);
   CHECK(basis.zeta_u == 0.01);
   CHECK(basis.zeta_d == 0.02);
   CHECK(basis.zeta_l == 0.03);
   CHECK(basis.yukawa_type == gm2calc::thdm::Yukawa_type::type_2);

   for (int i = 0, count = 1; i < 3; i++) {
      for (int k = 0; k < 3; k++, count++) {
         CHECK(std::real(basis.Xu(i,k)) == count);
         CHECK(std::real(basis.Xd(i,k)) == (count + 10));
         CHECK(std::real(basis.Xl(i,k)) == (count + 20));
         CHECK(std::imag(basis.Xu(i,k)) == (count + 30));
         CHECK(std::imag(basis.Xd(i,k)) == (count + 40));
         CHECK(std::imag(basis.Xl(i,k)) == (count + 50));
      }
   }
}
