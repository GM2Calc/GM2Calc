GM2Calc
=======

![](https://img.shields.io/github/v/release/GM2Calc/GM2Calc)
[![Build Status](https://github.com/GM2Calc/GM2Calc/workflows/test/badge.svg)](https://github.com/GM2Calc/GM2Calc/actions)

GM2Calc calculates the new physics contributions to the anomalous
magnetic moment of the muon a_mu = (g-2)/2 in the real MSSM and in the
THDM at the 1- and leading 2-loop level.

 * Homepage:                https://gm2calc.hepforge.org
 * Source code repository:  https://github.com/gm2calc
 * References:              Eur.Phys.J. C76 (2016) no.2, 62
   [[arxiv:1510.08071](https://arxiv.org/abs/1510.08071)]


Quick start
===========

Build GM2Calc as follows:

    mkdir -p build
    cd build
    cmake ..
    make

Run GM2Calc with input files at the command line:

    bin/gm2calc.x --gm2calc-input-file=../input/example.gm2
    bin/gm2calc.x --slha-input-file=../input/example.slha
    bin/gm2calc.x --thdm-input-file=../input/example.thdm


Contents
========

- [Requirements](#requirements)
- [Building GM2Calc](#building-gm2calc)
  * [Installation of dependencies](#installation-of-dependencies)
  * [Compilation of GM2Calc](#compilation-of-gm2calc)
- [Running GM2Calc](#running-gm2calc)
  * [From the command line](#from-the-command-line)
  * [From within Mathematica](#from-within-mathematica)
- [Input parameters](#input-parameters)
  * [MSSM: SLHA input parameters](#mssm-slha-input-parameters)
  * [MSSM: GM2Calc input parameters](#mssm-gm2calc-input-parameters)
  * [THDM: SLHA-like input parameters](#thdm-slha-like-input-parameters)
    + [Gauge basis](#gauge-basis)
    + [Mass basis](#mass-basis)
  * [Block `GM2CalcConfig`](#block-gm2calcconfig)
- [C/C++ interface](#cc-interface)
- [Mathematica interface](#mathematica-interface)
- [Source code documentation](#source-code-documentation)
- [References](#references)


Requirements
============

 * C++11 compatible compiler (g++ >= 4.6.3 or clang++ >= 3.4 or icpc >= 14.0.3)
 * C11 compatible compiler (gcc >= 4.6 or clang)
 * Boost (version 1.37.0 or higher) [http://www.boost.org]
 * Eigen 3 (version 3.1 or higher) [http://eigen.tuxfamily.org]

Optional:

 * Mathematica (version 9.0 or higher) [https://www.wolfram.com/mathematica/]


Building GM2Calc
================

Installation of dependencies
----------------------------

The dependencies listed above can be installed using the package
manager of your Linux distribution.  On Debian/Ubuntu, for example,
one may run:

    sudo apt-get install libeigen3-dev libboost-all-dev uuid-dev

Alternatively, the [Conan](https://conan.io/) package manager can be
used to install the dependencies:

    mkdir -p build
    cd build
    conan install ..

Compilation of GM2Calc
----------------------

To compile GM2Calc run:

    mkdir -p build
    cd build
    cmake ..
    make

The GM2Calc executable can then be found in `bin/gm2calc.x` and the
GM2Calc library can be found in the `lib/` directory.  The used
compiler and package paths can be passed as arguments to cmake.

**Example:**

    cmake \
       -DCMAKE_CXX_COMPILER=icpc \
       -DEIGEN3_INCLUDE_DIR=/opt/eigen3/eigen3 \
       -DBOOST_ROOT=/opt/boost \
       ..


Running GM2Calc
===============

From the command line
---------------------

GM2Calc can be run from the command line using an SLHA or SLHA-like
input. For the MSSM the input can be given in the SLHA format or in a
custom GM2Calc input format (similar to SLHA, but different definition
of input parameters in a mixed DR-bar/on-shell scheme). For the THDM
the input can be given in an SLHA-like format (compatible with
2HDMC). See `bin/gm2calc.x --help` for all options.

**Examples:** Running GM2Calc with an SLHA input file for the MSSM:

    bin/gm2calc.x --slha-input-file=../input/example.slha

or

    cat ../input/example.slha | bin/gm2calc.x --slha-input-file=-

**Example:** Running GM2Calc with a custom GM2Calc input file for the
MSSM (mixed DR-bar/on-shell scheme):

    bin/gm2calc.x --gm2calc-input-file=../input/example.gm2

or

    cat ../input/example.gm2 | bin/gm2calc.x --gm2calc-input-file=-

**Example:** Running GM2Calc with SLHA input file generated by
[SOFTSUSY](https://softsusy.hepforge.org/):

    ./softpoint.x leshouches < inOutFiles/lesHouchesInput | bin/gm2calc.x --slha-input-file=-

**Example:** Running GM2Calc with SLHA input file generated by
[SPheno](https://spheno.hepforge.org/):

    bin/SPheno input/LesHouches.in.mSUGRA
    bin/gm2calc.x --slha-input-file=SPheno.spc

**Example:** Running GM2Calc with an input file for the THDM:

    bin/gm2calc.x --thdminput-file=../input/example.thdm

or

    cat ../input/example.thdm | bin/gm2calc.x --thdm input-file=-


From within Mathematica
-----------------------

When Mathematica is installed, a MathLink executable for GM2Calc is
build, which allows one to run GM2Calc from within Mathematica.  The
MathLink executable can be found in `bin/gm2calc.mx`.  This executable
can be installed in Mathematica by calling

    Install["bin/gm2calc.mx"]

Afterwards, the GM2Calc Mathematica interface functions can be used.
See `examples/example-slha.m` (MSSM), `examples/example-gm2calc.m`
(MSSM) and `examples/example-thdm.m` (THDM) for examples with
different input parameters.

**Example:**

    math -run "<< ../examples/example-slha.m"
    math -run "<< ../examples/example-gm2calc.m"
    math -run "<< ../examples/example-thdm.m"

From within Python
------------------

When a valid python installation is detected, GM2Calc, will
attempt to create `example_*.py` files in `build/bin/`.  Either
Python2 or 3 can be used, and the python package `cppyy` is required.  
`cppyy` installation instructions can be found at:
[](https://cppyy.readthedocs.io/en/latest/installation.html)
See `examples/example_slha.py` (MSSM), `examples/example_gm2calc.py`
(MSSM) and `examples/example_thdm.py` (THDM) for examples with
different input parameters.

**Example:**

    python bin/example_slha.py
    python bin/example_gm2calc.py
    python bin/example_thdm.py


Input parameters
================

MSSM: SLHA input parameters
---------------------------

When GM2Calc is called with an
[SLHA-1](https://arxiv.org/abs/hep-ph/0311123) input file for the
MSSM, for example as

    bin/gm2calc.x --slha-input-file=../input/example.slha

then the input parameters are read from the following blocks:

 * `SMINPUTS`: α_s(MZ), Standard Model fermion masses, W and Z pole masses.

   *Example:*

       Block SMINPUTS
            3     0.1184           # alpha_s(MZ) SM MS-bar [2L]
            4     9.11876000E+01   # MZ(pole)              [1L]
            5     4.18000000E+00   # mb(mb) SM MS-bar      [2L]
            6     1.73340000E+02   # mtop(pole)            [2L]
            7     1.77700000E+00   # mtau(pole)            [2L]
            9     8.03850000E+01   # MW(pole)              [1L]
           13     0.1056583715     # mmuon(pole)           [1L]

 * `MASS`: pole masses of SUSY particles, W pole mass.  If the W pole
   mass is given in `MASS[24]`, then this value will be used instead
   of the W pole mass given in `SMINPUTS[9]`.

   *Example:*

       Block MASS                      # pole mass spectrum
               24     8.03773317e+01   # W                 [1L]
               36     1.50000000e+03   # A0                [2L]
          1000022     2.01611468e+02   # neutralino(1)     [1L]
          1000023     4.10040273e+02   # neutralino(2)     [1L]
          1000024     4.09989890e+02   # chargino(1)       [1L]
          1000025    -5.16529941e+02   # neutralino(3)     [1L]
          1000035     5.45628749e+02   # neutralino(4)     [1L]
          1000037     5.46057190e+02   # chargino(2)       [1L]
          1000013     5.25187016e+02   # smuon(1)          [1L]
          1000014     5.18860573e+02   # muon sneutrino    [1L]
          2000013     5.05095249e+02   # smuon(2)          [1L]


 * `HMIX`: μ parameter (initial guess), tan(β), renormalization scale Q

   *Example:*

       Block HMIX Q= 1.00000000e+03   # renormalization scale
            1     4.89499929e+02      # mu(Q) MSSM DR-bar   [initial guess]
            2     3.93371545e+01      # tan(beta)(Q) MSSM DR-bar Feynman gauge [1L]

 * `AU`, `AD`, `AE`: soft-breaking trilinear couplings at the
   renormalization scale Q, which has been read from the `HMIX` block

   *Example:*

       Block AU Q= 1.00000000e+03
         3  3     1.57871614e-05      # At(Q)              [2L]
       Block AD Q= 1.00000000e+03
         3  3     8.99561673e-06      # Ab(Q)              [2L]
       Block AE Q= 1.00000000e+03
         2  2     2.84230475e-06      # Amu(Q) MSSM DR-bar [1L]
         3  3     3.02719242e-06      # Atau(Q)            [2L]

 * `MSOFT`: soft-breaking gaugino masses, soft-breaking sfermion
   masses at the renormalization scale Q, which has been read from the
   `HMIX` block

   *Example:*

       Block MSOFT Q= 1.00000000e+03
            1     2.00000000e+02      # M_1(Q)             [initial guess]
            2     4.00000000e+02      # M_2(Q)             [initial guess]
            3     2.00000000e+03      # M_3(Q)             [2L]
           31     5.00000000e+02      # meL(Q)             [2L]
           32     5.00000000e+02      # mmuL(Q)            [irrelevant]
           33     3.00000000e+03      # mtauL(Q)           [2L]
           34     4.99999999e+02      # meR(Q)             [2L]
           35     4.99999999e+02      # mmuR(Q)            [initial guess]
           36     3.00000000e+03      # mtauR(Q)           [2L]
           41     7.00000000e+03      # mqL1(Q)            [2L]
           42     7.00000000e+03      # mqL2(Q)            [2L]
           43     6.99999999e+03      # mqL3(Q)            [2L]
           44     7.00000000e+03      # muR(Q)             [2L]
           45     7.00000000e+03      # mcR(Q)             [2L]
           46     6.99999999e+03      # mtR(Q)             [2L]
           47     7.00000000e+03      # mdR(Q)             [2L]
           48     7.00000000e+03      # msR(Q)             [2L]
           49     7.00000000e+03      # mbR(Q)             [2L]

 * `GM2CalcInput`: α_em(MZ), α_em(0) in the Thomson limit

   *Example:*

       Block GM2CalcInput
            1     0.00775520       # alpha(MZ)             [1L]
            2     0.00729735       # alpha(0)              [2L]

See `input/example.slha` for an example SLHA input file.

---
**Important note:**

We strongly *encourage* passing SLHA input in
[SLHA-1](https://arxiv.org/abs/hep-ph/0311123) convention to GM2Calc.
We strongly *discourage* passing
[SLHA-2](https://arxiv.org/abs/0801.0045) input to GM2Calc. The reason
is, that GM2Calc interprets the PDG numbers 1000013, 2000013 and
1000014 in the `MASS` block as the two smuon and the muon sneutrino
pole masses, respectively. In an SLHA-2-compliant input, however, it
is not ensured, that these PDG numbers correspond to the two smuon and
the muon sneutrino pole masses, because of possibly allowed slepton
flavour violation.

---

MSSM: GM2Calc input parameters
------------------------------

When GM2Calc is called with an input file in the custom GM2Calc
format (in a mixed DR-bar/on-shell scheme), for example as

    bin/gm2calc.x --gm2calc-input-file=../input/example.gm2

then the input parameters are read from the following blocks:

 * `SMINPUTS`: α_s(MZ), Standard Model fermion masses, W and Z pole masses.

   *Example:*

       Block SMINPUTS
            3     0.1184           # alpha_s(MZ) SM MS-bar [2L]
            4     9.11876000E+01   # MZ(pole)              [1L]
            5     4.18000000E+00   # mb(mb) SM MS-bar      [2L]
            6     1.73340000E+02   # mtop(pole)            [2L]
            7     1.77700000E+00   # mtau(pole)            [2L]
            9     8.03850000E+01   # MW(pole)              [1L]
           13     0.1056583715     # mmuon(pole)           [1L]

 * `GM2CalcInput`: renormalization scale Q, α_em(MZ), α_em(0) in the
   Thomson limit, tan(β) (DR-bar scheme), μ parameter (on-shell),
   soft-breaking gaugino masses (M1 and M2 in on-shell scheme), CP-odd
   Higgs pole mass MA, soft-breaking sfermion masses (msl(2,2) and
   mse(2,2) in on-shell scheme), soft-breaking trilinear couplings
   (Ae(2,2) in DR-bar scheme)

   *Example:*

       Block GM2CalcInput
            0     8.66360379E+02   # ren. scale Q          [2L]
            1     0.00775531       # alpha(MZ)             [1L]
            2     0.00729735       # alpha(0)              [2L]
            3     10               # tan(beta) DR-bar at Q [1L]
            4     619.858          # mu parameter on-shell [1L]
            5     211.722          # M1 on-shell           [1L]
            6     401.057          # M2 on-shell           [1L]
            7     1.10300877E+03   # M3                    [2L]
            8     707.025          # MA(pole)              [2L]
            9     3.51653258E+02   # msl(1,1)              [2L]
           10     356.09           # msl(2,2) on-shell     [1L]
           11     3.50674223E+02   # msl(3,3)              [2L]
           12     2.21215037E+02   # mse(1,1)              [2L]
           13     225.076          # mse(2,2) on-shell     [1L]
           14     2.18022142E+02   # mse(3,3)              [2L]
           15     1.00711403E+03   # msq(1,1)              [2L]
           16     1.00711149E+03   # msq(2,2)              [2L]
           17     9.29083096E+02   # msq(3,3)              [2L]
           18     9.69369660E+02   # msu(1,1)              [2L]
           19     9.69366965E+02   # msu(2,2)              [2L]
           20     7.99712943E+02   # msu(3,3)              [2L]
           21     9.64756473E+02   # msd(1,1)              [2L]
           22     9.64753818E+02   # msd(2,2)              [2L]
           23     9.60016201E+02   # msd(3,3)              [2L]
           24     0                # Ae(1,1)               [irrelevant]
           25    -2.93720212E+02   # Ae(2,2) DR-bar        [1L]
           26    -2.92154796E+02   # Ae(3,3)               [2L]
           27     0                # Ad(1,1)               [irrelevant]
           28     0                # Ad(2,2)               [irrelevant]
           29    -1.28330100E+03   # Ad(3,3)               [2L]
           30     0                # Au(1,1)               [irrelevant]
           31     0                # Au(2,2)               [irrelevant]
           32    -8.70714986E+02   # Au(3,3)               [2L]

See `input/example.gm2` for an example input file with custom GM2Calc
input parameters.

THDM: SLHA-like input parameters
--------------------------------

When GM2Calc is called with an SLHA-like input file for the THDM, for
example as

    bin/gm2calc.x --thdm-input-file=../input/example.thdm

then the input parameters are read from different blocks, depending on
the input basis.

### Gauge basis

The THDM-specific input parameters in the gauge basis are read from
the following blocks:

 * `MINPAR`: tan(β), λ\_{1,...,7}, m\_{12}^2, ζ\_{u,d,l} and the
   Yukawa type (1 = type I, 2 = type II, 3 = type X, 4 = type Y, 5 =
   aligned THDM, 6 = general)

   Note: The parameters ζ\_{u,d,l} are only used if the Yukawa type is
   set to the aligned THDM (Yukawa type = 5).

   *Example:*

       Block MINPAR                  # model parameters in gauge basis
           3        3                # tan(beta)
          11        0.7              # lambda_1
          12        0.6              # lambda_2
          13        0.5              # lambda_3
          14        0.4              # lambda_4
          15        0.3              # lambda_5
          16        0.2              # lambda_6
          17        0.1              # lambda_7
          18        40000            # m_{12}^2
          21        0                # zeta_u
          22        0                # zeta_d
          23        0                # zeta_l
          24        2                # Yukawa type (1, 2, 3, 4, 5 = aligned, 6 = general)

 * `GM2CalcTHDMDeltauInput`, `GM2CalcTHDMDeltadInput`, `GM2CalcTHDMDeltalInput`:
   The real parts of the matrices Delta_u, Delta_d and Delta_l

 * `GM2CalcTHDMPiuInput`, `GM2CalcTHDMPidInput`, `GM2CalcTHDMPilInput`:
   The real parts of the matrices Pi_u, Pi_d and Pi_l

Note: The parameters Delta\_{u,d,l} are only used if the Yukawa type is
not set to the general THDM (Yukawa type = 1,...,5).

Note: The parameters Pi\_{u,d,l} are only used if the Yukawa type is
set to the general THDM (Yukawa type = 6).

### Mass basis

The THDM-specific input parameters in the mass basis are read from the
following blocks:

 * `MINPAR`: tan(β), λ\_{6,7}, m\_{12}^2, sin(β-α), ζ\_{u,d,l} and the
   Yukawa type (1 = type I, 2 = type II, 3 = type X, 4 = type Y, 5 =
   aligned THDM, 6 = general)

   Note: The parameters ζ\_{u,d,l} are only used if the Yukawa type is
   set to the aligned THDM (Yukawa type = 5).

   *Example:*

       Block MINPAR                  # model parameters
           3        3                # tan(beta)
          16        0.2              # lambda_6
          17        0.1              # lambda_7
          18        40000            # m_{12}^2
          20        0.999            # sin(beta - alpha)
          21        0                # zeta_u
          22        0                # zeta_d
          23        0                # zeta_l
          24        2                # Yukawa type (1, 2, 3, 4, 5 = aligned, 6 = general)

 * `MASS`: CP-even Higgs boson masses m_h, m_H, the CP-odd Higgs boson
   mass m_A and the charged Higgs boson mass m_{H^+}

   *Example:*

       Block MASS                    # Higgs masses
          25        100              # mh, lightest CP-even Higgs
          35        400              # mH, heaviest CP-even Higgs
          36        420              # mA, CP-odd Higgs
          37        440              # mH+, charged Higgs

 * `GM2CalcTHDMDeltauInput`, `GM2CalcTHDMDeltadInput`, `GM2CalcTHDMDeltalInput`:
   The real parts of the matrices Delta_u, Delta_d and Delta_l

 * `GM2CalcTHDMPiuInput`, `GM2CalcTHDMPidInput`, `GM2CalcTHDMPilInput`:
   The real parts of the matrices Pi_u, Pi_d and Pi_l

Note: The parameters Delta\_{u,d,l} are only used if the Yukawa type is
not set to the general THDM (Yukawa type = 1,...,5).

Note: The parameters Pi\_{u,d,l} are only used if the Yukawa type is
set to the general THDM (Yukawa type = 6).

See `input/example.thdm` for an example input file.


Block `GM2CalcConfig`
---------------------

When running GM2Calc from the command line with SLHA or SLHA-like
input, the input may contain the `GM2CalcConfig` configuration block
to customize the calculation and the output.  The `GM2CalcConfig`
block entries are summarized in the following table and are described
below:

| Entry    | Description           | Possible values | Defaul value                              |
|----------|-----------------------|-----------------|-------------------------------------------|
| 0        | output format         | `0`, ..., `4`   | `4` for SLHA input, `1` for GM2Calc input |
| 1        | loop order            | `0`, `1`, `2`   | `2`                                       |
| 2        | tan(beta) resummation | `0`, `1`        | `1`                                       |
| 3        | force output          | `0`, `1`        | `0`                                       |
| 4        | verboe output         | `0`, `1`        | `0`                                       |
| 5        | estimate uncertainty  | `0`, `1`        | `0`                                       |
| 6        | running parameters    | `0`, `1`        | `1`                                       |

Description:

 * `GM2CalcConfig[0]`: defines the output format (`0` ... `4`)

   `0`: minimal output (a single number)  
   `1`: detailed (a detailed output of the various contributions)  
   `2`: write the value of a_mu into `LOWEN` block, entry `6`  
   `3`: write the value of a_mu into `SPhenoLowEnergy` block, entry `21`  
   `4`: write the value of a_mu into `GM2CalcOutput` block, entry `0`

 * `GM2CalcConfig[1]`: loop order of the calculation (`0`, `1` or `2`).
   We recommend to use `2`.

 * `GM2CalcConfig[2]`: disable/enable tan(beta) resummation (`0` or `1`).
   We recommend to use `1`.

 * `GM2CalcConfig[3]`: force output even if physical problem has
   occured (`0` or `1`).  **WARNING**: The result might not be trusted
   if a problem has occured!  We recommend to use `0`.

 * `GM2CalcConfig[4]`: disable/enable verbose output (`0` or `1`).
   We recommend to use `0` by default, and `1` only for debugging.

 * `GM2CalcConfig[5]`: disable/enable uncertainty estimation (`0` or `1`).
   Depending on the chosen output format (see `GM2CalcConfig[0]`) the
   uncertainty is written

   * as a single number to stdout   in case of minimal output,
   * to the first line              in case of detailed output,
   * to `GM2CalcOutput[1]`          otherwise.

 * `GM2CalcConfig[6]`: disable/enable use of running parameters (`0` or `1`)
   in the fermionic 2-loop contributions in the THDM.

We recommend to use the following configuration block in an SLHA input
file:

    Block GM2CalcConfig
         0     4     # output format (0 = minimal, 1 = detailed,
                     #  2 = NMSSMTools, 3 = SPheno, 4 = GM2Calc)
         1     2     # loop order (0, 1 or 2)
         2     1     # disable/enable tan(beta) resummation (0 or 1)
         3     0     # force output (0 or 1)
         4     0     # verbose output (0 or 1)
         5     1     # calculate uncertainty (0 or 1)
         6     1     # running parameters (0 or 1)


C/C++ interface
===============

GM2Calc provides a C and a C++ interface.  To use the routines of
GM2Calc in a C++ program and perform an MSSM calculation, the
following C++ header files have to be included:

    include/gm2calc/gm2_1loop.hpp
    include/gm2calc/gm2_2loop.hpp
    include/gm2calc/gm2_uncertainty.hpp
    include/gm2calc/gm2_error.hpp
    include/gm2calc/MSSMNoFV_onshell.hpp

For the calculation in the THDM, the following C++ header files have
to be included:

    include/gm2calc/gm2_1loop.hpp
    include/gm2calc/gm2_2loop.hpp
    include/gm2calc/gm2_uncertainty.hpp
    include/gm2calc/gm2_error.hpp
    include/gm2calc/THDM.hpp

To use the routines of GM2Calc in a C program to perform an MSSM
calculation, the following C header files have to be included:

    include/gm2calc/gm2_1loop.h
    include/gm2calc/gm2_2loop.h
    include/gm2calc/gm2_uncertainty.h
    include/gm2calc/MSSMNoFV_onshell.h

For the calculation in the THDM, the following C header files must be
included:

    include/gm2calc/gm2_1loop.h
    include/gm2calc/gm2_2loop.h
    include/gm2calc/gm2_uncertainty.h
    include/gm2calc/THDM.h
    include/gm2calc/SM.h

Please refer to the content of these header files for a precise
definition of all interface functions.  The C/C++ example programs in
the `examples/` directory serve as an illustration of the interface
routines.


Mathematica interface
=====================

After the GM2Calc MathLink executable has been loaded by
e.g. `Install["bin/gm2calc.mx"]`, the following functions are
available:

| Function                | Description                                                           |
| ----------------------- | --------------------------------------------------------------------- |
| GM2CalcSetFlags         | Sets configuration flags for the calculation                          |
| GM2CalcGetFlags         | Returns currently set configuration flags                             |
| GM2CalcSetSMParameters  | Sets Standard Model parameters                                        |
| GM2CalcGetSMParameters  | Returns currently set Standard Model parameters                       |
| GM2CalcAmuSLHAScheme    | Calculates `a_mu`, in the MSSM with SLHA input parameters             |
| GM2CalcAmuGM2CalcScheme | Calculates `a_mu`, in the MSSM with GM2Calc-specific input parameters |
| GM2CalcAmuTHDMGaugeBasis| Calculates `a_mu`, in the THDM with gauge basis input parameters      |
| GM2CalcAmuTHDMMassBasis | Calculates `a_mu`, in the THDM with mass basis input parameters       |

See the example Mathematica scripts `examples/example-slha.m`,
`examples/example-gm2calc.m` and `examples/example-thdm.m`.


Python interface
================

After building GM2Calc with shared libraries, it is possible to
load gm2calc functions into python using the C-Python interface
provided by cppyy.  

To use the routines of GM2Calc in a python script, the following 
folders and C++ header files have to be included:

    `GM2Calc include folder`
    `Eigen3 include folder`
    include/gm2calc/gm2_1loop.hpp
    include/gm2calc/gm2_2loop.hpp
    include/gm2calc/gm2_uncertainty.hpp
    include/gm2calc/gm2_error.hpp

To perform an MSSM calculation, the following C++ header files 
have to be included:

    include/gm2calc/MSSMNoFV_onshell.hpp

For the calculation in the THDM, the following C++ header files 
have to be included:

    include/gm2calc/THDM.hpp

And finally for any GM2Calc model the following library must be
loaded:

    cppyy.load_library("libgm2calc")

Then one can import the C++ functions into python using the command
`from cppy.gbl import gm2calc`.  Then the following functions are 
available:

| Function                | Description                                                           |
| ----------------------- | --------------------------------------------------------------------- |
| gm2calc.SM              | Constructs an SM model                                                |
| gm2calc.MSSM_onshell    | Constructs an MSSM model                                              |
| gm2calc.THDM            | Constructs a THDM model                                               |
| gm2calc.calculate_uncertainty_amu_0loop | Calculates the uncertainty in contributions to `a_mu` if working at tree-level |
| gm2calc.calculate_uncertainty_amu_1loop | Calculates the uncertainty in contributions to `a_mu` if working at 1 level    |
| gm2calc.calculate_uncertainty_amu_2loop | Calculates the uncertainty in contributions to `a_mu` if working at 2 level    |

See the example Python scripts `examples/example_slha.py`,
`examples/example_gm2calc.py` and `examples/example_thdm.py`.


Source code documentation
=========================

GM2Calc uses Doxygen to document the source code.  The complete source
code documentation can be generated by running

    make doc

Afterwards, the generated HTML pages can be found in `doc/html/` and
can be openend with your favourite web browser, e.g.

    firefox doc/html/index.html


References
==========

GM2Calc has been published in [`Athron:2015rva`]

- Peter Athron, Markus Bach, Helvecio G. Fargnoli, Christoph
  Gnendiger, Robert Greifenhagen, Jae-hyeon Park, Sebastian Paßehr,
  Dominik Stöckinger, Hyejung Stöckinger-Kim, Alexander Voigt:
  *GM2Calc: Precise MSSM prediction for (g-2) of the muon* [<a
  href="http://arxiv.org/abs/1510.08071">Eur.Phys.J. C76 (2016) no.2,
  62</a>]

- Peter Athron, Csaba Balazs, Adriano Cherchiglia, Douglas Jacob,
  Dominik Stöckinger, Hyejung Stöckinger-Kim, Alexander Voigt:
  *[TODO]*

The expressions implemented in GM2Calc for the MSSM have been taken from
[<a href="http://arxiv.org/abs/1003.5820">Phys.Rev. D81 (2010) 093004</a>,
 <a href="http://arxiv.org/abs/1309.0980">Phys.Lett. B726 (2013) 717-724</a>,
 <a href="http://arxiv.org/abs/1311.1775">JHEP 1402 (2014) 070</a>,
 <a href="http://arxiv.org/abs/1504.05500">JHEP 1510 (2015) 026</a>]

The expressions implemented in GM2Calc for the THDM have been taken from
[<a href="http://arxiv.org/abs/1607.06292">JHEP 01 (2017) 007</a>]
