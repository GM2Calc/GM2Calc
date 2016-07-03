DIR          := src
MODNAME      := gm2calc

# source files for library
LIBsrc_SRC := \
		$(DIR)/dilog.cpp \
		$(DIR)/ffunctions.cpp \
		$(DIR)/gm2_1loop.cpp \
		$(DIR)/gm2_1loop_c.cpp \
		$(DIR)/gm2_2loop.cpp \
		$(DIR)/gm2_2loop_c.cpp \
		$(DIR)/gm2_mb.cpp \
		$(DIR)/gm2_uncertainty.cpp \
		$(DIR)/gm2_uncertainty_c.cpp \
		$(DIR)/gm2_slha_io.cpp \
		$(DIR)/MSSMNoFV_onshell.cpp \
		$(DIR)/MSSMNoFV_onshell_c.cpp \
		$(DIR)/MSSMNoFV_onshell_mass_eigenstates.cpp \
		$(DIR)/MSSMNoFV_onshell_physical.cpp \
		$(DIR)/MSSMNoFV_onshell_problems.cpp \
		$(DIR)/MSSMNoFV_onshell_soft_parameters.cpp \
		$(DIR)/MSSMNoFV_onshell_susy_parameters.cpp

# source files with main()
EXEsrc_SRC := \
		$(DIR)/gm2calc.cpp

MATHLINK_SRC := $(DIR)/gm2calc.tm

LIBsrc_OBJ := $(LIBsrc_SRC:.cpp=.o)

EXEsrc_OBJ := $(EXEsrc_SRC:.cpp=.o)

LIBsrc_DEP := $(LIBsrc_OBJ:.o=.d)

EXEsrc_DEP := $(EXEsrc_OBJ:.o=.d)

EXEsrc_EXE := $(patsubst $(DIR)/%.o, $(BINDIR)/%.x, $(EXEsrc_OBJ))

MATHLINK_EXE := $(patsubst $(DIR)/%.tm, $(BINDIR)/%.mx, $(MATHLINK_SRC))

LIBsrc     := $(DIR)/lib$(MODNAME)$(LIBEXT)

SHAREDLIBsrc := $(DIR)/lib$(MODNAME)$(SHAREDLIBEXT)

.PHONY: mathlink

clean::
	-rm -f $(LIBsrc_DEP) $(EXEsrc_DEP)
	-rm -f $(LIBsrc_OBJ) $(EXEsrc_OBJ)
	-rm -f $(LIBsrc) $(SHAREDLIBsrc)
	-rm -f $(EXEsrc_EXE)
	-rm -f $(MATHLINK_EXE)

mathlink: $(BINDIR)/gm2calc.mx

$(LIBsrc_DEP) $(EXEsrc_DEP) $(LIBsrc_OBJ) $(EXEsrc_OBJ): \
	override CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS)

$(LIBsrc): $(LIBsrc_OBJ)
	$(MAKELIB) $@ $^

$(SHAREDLIBsrc): $(LIBsrc_OBJ)
	$(MAKESHAREDLIB) -o $@ $^

$(BINDIR)/%.x: $(DIR)/%.o $(LIBsrc)
	$(CXX) -o $@ $^ $(LDLIBS)

$(BINDIR)/%.mx: $(DIR)/%.tm $(LIBsrc) | $(FCC) $(FXX)
	CXX="$(FXX)" REALCXX="$(CXX) $(CXXFLAGS)" \
	CC="$(FCC)" REALCC="$(CC) $(CFLAGS)" MATH="$(MATH)" \
	$(MCC) -o $@ $(CPPFLAGS) $(CFLAGS) $^ $(LDLIBS)

ALLDEP += $(LIBsrc_DEP) $(EXEsrc_DEP)
ALLLIB += $(LIBsrc) $(SHAREDLIBsrc)
ALLEXE += $(EXEsrc_EXE)
