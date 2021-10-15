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
Block GM2CalcTHDMDeltauInput
    1 1   1                    # Re(Delta_u(1,1))
    1 2   2                    # Re(Delta_u(1,2))
    1 3   3                    # Re(Delta_u(1,3))
    2 1   4                    # Re(Delta_u(2,1))
    2 2   5                    # Re(Delta_u(2,2))
    2 3   6                    # Re(Delta_u(2,3))
    3 1   7                    # Re(Delta_u(3,1))
    3 2   8                    # Re(Delta_u(3,2))
    3 3   9                    # Re(Delta_u(3,3))
Block GM2CalcTHDMDeltadInput
    1 1   11                   # Re(Delta_d(1,1))
    1 2   12                   # Re(Delta_d(1,2))
    1 3   13                   # Re(Delta_d(1,3))
    2 1   14                   # Re(Delta_d(2,1))
    2 2   15                   # Re(Delta_d(2,2))
    2 3   16                   # Re(Delta_d(2,3))
    3 1   17                   # Re(Delta_d(3,1))
    3 2   18                   # Re(Delta_d(3,2))
    3 3   19                   # Re(Delta_d(3,3))
Block GM2CalcTHDMDeltalInput
    1 1   21                   # Re(Delta_l(1,1))
    1 2   22                   # Re(Delta_l(1,2))
    1 3   23                   # Re(Delta_l(1,3))
    2 1   24                   # Re(Delta_l(2,1))
    2 2   25                   # Re(Delta_l(2,2))
    2 3   26                   # Re(Delta_l(2,3))
    3 1   27                   # Re(Delta_l(3,1))
    3 2   28                   # Re(Delta_l(3,2))
    3 3   29                   # Re(Delta_l(3,3))
Block GM2CalcTHDMPiuInput
    1 1   31                   # Re(Pi_u(1,1))
    1 2   32                   # Re(Pi_u(1,2))
    1 3   33                   # Re(Pi_u(1,3))
    2 1   34                   # Re(Pi_u(2,1))
    2 2   35                   # Re(Pi_u(2,2))
    2 3   36                   # Re(Pi_u(2,3))
    3 1   37                   # Re(Pi_u(3,1))
    3 2   38                   # Re(Pi_u(3,2))
    3 3   39                   # Re(Pi_u(3,3))
Block GM2CalcTHDMPidInput
    1 1   41                   # Re(Pi_d(1,1))
    1 2   42                   # Re(Pi_d(1,2))
    1 3   43                   # Re(Pi_d(1,3))
    2 1   44                   # Re(Pi_d(2,1))
    2 2   45                   # Re(Pi_d(2,2))
    2 3   46                   # Re(Pi_d(2,3))
    3 1   47                   # Re(Pi_d(3,1))
    3 2   48                   # Re(Pi_d(3,2))
    3 3   49                   # Re(Pi_d(3,3))
Block GM2CalcTHDMPilInput
    1 1   51                   # Re(Pi_l(1,1))
    1 2   52                   # Re(Pi_l(1,2))
    1 3   53                   # Re(Pi_l(1,3))
    2 1   54                   # Re(Pi_l(2,1))
    2 2   55                   # Re(Pi_l(2,2))
    2 3   56                   # Re(Pi_l(2,3))
    3 1   57                   # Re(Pi_l(3,1))
    3 2   58                   # Re(Pi_l(3,2))
    3 3   59                   # Re(Pi_l(3,3))
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
         CHECK(basis.Delta_u(i,k) == count);
         CHECK(basis.Delta_d(i,k) == (count + 10));
         CHECK(basis.Delta_l(i,k) == (count + 20));
         CHECK(basis.Pi_u(i,k) == (count + 30));
         CHECK(basis.Pi_d(i,k) == (count + 40));
         CHECK(basis.Pi_l(i,k) == (count + 50));
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
Block GM2CalcTHDMDeltauInput
    1 1   1                    # Re(Delta_u(1,1))
    1 2   2                    # Re(Delta_u(1,2))
    1 3   3                    # Re(Delta_u(1,3))
    2 1   4                    # Re(Delta_u(2,1))
    2 2   5                    # Re(Delta_u(2,2))
    2 3   6                    # Re(Delta_u(2,3))
    3 1   7                    # Re(Delta_u(3,1))
    3 2   8                    # Re(Delta_u(3,2))
    3 3   9                    # Re(Delta_u(3,3))
Block GM2CalcTHDMDeltadInput
    1 1   11                   # Re(Delta_d(1,1))
    1 2   12                   # Re(Delta_d(1,2))
    1 3   13                   # Re(Delta_d(1,3))
    2 1   14                   # Re(Delta_d(2,1))
    2 2   15                   # Re(Delta_d(2,2))
    2 3   16                   # Re(Delta_d(2,3))
    3 1   17                   # Re(Delta_d(3,1))
    3 2   18                   # Re(Delta_d(3,2))
    3 3   19                   # Re(Delta_d(3,3))
Block GM2CalcTHDMDeltalInput
    1 1   21                   # Re(Delta_l(1,1))
    1 2   22                   # Re(Delta_l(1,2))
    1 3   23                   # Re(Delta_l(1,3))
    2 1   24                   # Re(Delta_l(2,1))
    2 2   25                   # Re(Delta_l(2,2))
    2 3   26                   # Re(Delta_l(2,3))
    3 1   27                   # Re(Delta_l(3,1))
    3 2   28                   # Re(Delta_l(3,2))
    3 3   29                   # Re(Delta_l(3,3))
Block GM2CalcTHDMPiuInput
    1 1   31                   # Re(Pi_u(1,1))
    1 2   32                   # Re(Pi_u(1,2))
    1 3   33                   # Re(Pi_u(1,3))
    2 1   34                   # Re(Pi_u(2,1))
    2 2   35                   # Re(Pi_u(2,2))
    2 3   36                   # Re(Pi_u(2,3))
    3 1   37                   # Re(Pi_u(3,1))
    3 2   38                   # Re(Pi_u(3,2))
    3 3   39                   # Re(Pi_u(3,3))
Block GM2CalcTHDMPidInput
    1 1   41                   # Re(Pi_d(1,1))
    1 2   42                   # Re(Pi_d(1,2))
    1 3   43                   # Re(Pi_d(1,3))
    2 1   44                   # Re(Pi_d(2,1))
    2 2   45                   # Re(Pi_d(2,2))
    2 3   46                   # Re(Pi_d(2,3))
    3 1   47                   # Re(Pi_d(3,1))
    3 2   48                   # Re(Pi_d(3,2))
    3 3   49                   # Re(Pi_d(3,3))
Block GM2CalcTHDMPilInput
    1 1   51                   # Re(Pi_l(1,1))
    1 2   52                   # Re(Pi_l(1,2))
    1 3   53                   # Re(Pi_l(1,3))
    2 1   54                   # Re(Pi_l(2,1))
    2 2   55                   # Re(Pi_l(2,2))
    2 3   56                   # Re(Pi_l(2,3))
    3 1   57                   # Re(Pi_l(3,1))
    3 2   58                   # Re(Pi_l(3,2))
    3 3   59                   # Re(Pi_l(3,3))
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
   CHECK(basis.lambda_6 == 0.2);
   CHECK(basis.lambda_7 == 0.1);
   CHECK(basis.sin_beta_minus_alpha == 0.995);
   CHECK(basis.m122 == 40000.0);
   CHECK(basis.zeta_u == 0.01);
   CHECK(basis.zeta_d == 0.02);
   CHECK(basis.zeta_l == 0.03);
   CHECK(basis.yukawa_type == gm2calc::thdm::Yukawa_type::type_2);

   for (int i = 0, count = 1; i < 3; i++) {
      for (int k = 0; k < 3; k++, count++) {
         CHECK(basis.Delta_u(i,k) == count);
         CHECK(basis.Delta_d(i,k) == (count + 10));
         CHECK(basis.Delta_l(i,k) == (count + 20));
         CHECK(basis.Pi_u(i,k) == (count + 30));
         CHECK(basis.Pi_d(i,k) == (count + 40));
         CHECK(basis.Pi_l(i,k) == (count + 50));
      }
   }
}
