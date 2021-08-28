#!/usr/bin/env python

from __future__ import print_function

import cppyy

models = ["MSSMNoFV","THDM","SM"]

gm2calcpath = "@PROJECT_SOURCE_DIR@"
# Add gm2calc headers to path
includepath = gm2calcpath + "/include"
sourcepath = gm2calcpath + "/src"
cppyy.add_include_path(includepath)
cppyy.add_include_path(sourcepath)
for modelname in models:
	cppyy.add_include_path(sourcepath+"/"+modelname)
    
# Add Eigen3/Core directory to path
eigen3path = "@EIGEN3_INCLUDE_DIR@"
cppyy.add_include_path(eigen3path)

cppyy.include("gm2calc/gm2_1loop.hpp")
cppyy.include("gm2calc/gm2_2loop.hpp")
cppyy.include("gm2calc/gm2_uncertainty.hpp")
cppyy.include("gm2calc/gm2_error.hpp")
cppyy.include("gm2calc/THDM.hpp")

# Load library containing gm2calc definitions
librarypath = "@CMAKE_LIBRARY_OUTPUT_DIRECTORY@"
cppyy.add_library_path(librarypath)
cppyy.load_library("libgm2calc")

# Load data types
from cppyy.gbl import std
from cppyy.gbl import Eigen
# Shorten class call name
Matrix3cd = Eigen.Matrix[std.complex['double'],3,3]

from cppyy.gbl import gm2calc
from cppyy.gbl.gm2calc import SM
from cppyy.gbl.gm2calc import THDM
from cppyy.gbl.gm2calc import Error
basis = gm2calc.thdm.Mass_basis()

basis.yukawa_type = gm2calc.thdm.Yukawa_type.type_2
basis.mh = 125.
basis.mH = 400.
basis.mA = 420.
basis.mHp = 440.
basis.sin_beta_minus_alpha = 0.999
basis.lambda_6 = 0.
basis.lambda_7 = 0.
basis.tan_beta = 3.
basis.m122 = 40000.
basis.zeta_u = 0.
basis.zeta_d = 0.
basis.zeta_l = 0.
basis.Xu = Matrix3cd().setZero()
basis.Xd = Matrix3cd().setZero()
basis.Xl = Matrix3cd().setZero()

sm = gm2calc.SM()
sm.set_alpha_em_mz(1.0/128.94579)
sm.set_mu(2,173.34)
sm.set_mu(1,1.28)
sm.set_md(2,4.18)
sm.set_ml(2,1.77684)

config = gm2calc.thdm.Config()
config.force_output = False;
config.running_couplings = True;

try:
    model = gm2calc.THDM(basis,sm,config)
    amu = gm2calc.calculate_amu_1loop(model) + gm2calc.calculate_amu_2loop(model)
    delta_amu = gm2calc.calculate_uncertainty_amu_2loop(model)
    print("amu=",amu,"+-",delta_amu)
except gm2calc.Error as e:
    print(e.what)

