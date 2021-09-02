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
cppyy.include("gm2calc/MSSMNoFV_onshell.hpp")

# Load library containing gm2calc definitions
librarypath = "@CMAKE_LIBRARY_OUTPUT_DIRECTORY@"
cppyy.add_library_path(librarypath)
cppyy.load_library("libgm2calc")	
	
# Load data types
from cppyy.gbl import gm2calc
from cppyy.gbl.gm2calc import Error	

def setup():
    # load data types
    from cppyy.gbl import std
    from cppyy.gbl import Eigen
    # Shorten class call name
    from cppyy.gbl.Eigen import Matrix3d
    
    model = gm2calc.MSSMNoFV_onshell()
    
    Pi = 3.141592653589793
    # Outer Matrix3cd is to convert type from CwiseNullaryOp to Matrix
    UnitMatrix = Matrix3d(Matrix3d().Identity())
    # __mul__ is not defined for Matrix type
    UnitMatrix.__imul__(500*500)
    
    # fill SM parameters
    model.set_alpha_MZ(0.0077552)               # 1L
    model.set_alpha_thompson(0.00729735)        # 2L
    model.set_g3(std.sqrt(4. * Pi * 0.1184))    # 2L
    model.get_physical().MFt   = 173.34         # 2L
    model.get_physical().MFb   = 4.18           # 2L, mb(mb) MS-bar
    model.get_physical().MFm   = 0.1056583715   # 1L
    model.get_physical().MFtau = 1.777          # 2L
    model.get_physical().MVWm  = 80.385         # 1L
    model.get_physical().MVZ   = 91.1876        # 1L

    # fill DR-bar parameters
    model.set_TB(10)                            # 1L
    model.set_Ae(1,1,0)                         # 1L
    
    # fill on-shell parameters
    model.set_Mu(350)                      # 1L
    model.set_MassB(150)                   # 1L
    model.set_MassWB(300)                  # 1L
    model.set_MassG(1000)                  # 2L
    model.set_mq2(UnitMatrix)              # 2L
    model.set_ml2(UnitMatrix)              # 1L(smuon)/2L
    model.set_md2(UnitMatrix)              # 2L
    model.set_mu2(UnitMatrix)              # 2L
    model.set_me2(UnitMatrix)              # 1L(smuon)/2L
    model.set_Au(2,2,0)                    # 2L
    model.set_Ad(2,2,0)                    # 2L
    model.set_Ae(2,2,0)                    # 2L
    model.set_MA0(1500)                    # 2L
    model.set_scale(454.7)                 # 2L

    # calculate mass spectrum
    model.calculate_masses()
    
    return model


try:
    model = gm2calc.MSSMNoFV_onshell(setup())
    amu = gm2calc.calculate_amu_1loop(model) + gm2calc.calculate_amu_2loop(model)
    delta_amu = gm2calc.calculate_uncertainty_amu_2loop(model)
    print("amu =",amu,"+-",delta_amu)
except gm2calc.Error as e:
    print(e.what)
