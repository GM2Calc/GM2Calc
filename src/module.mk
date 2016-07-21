DIR          := src
MODNAME      := gm2calc

# source files with main()
EXEsrc_SRC := $(DIR)/gm2calc.cpp

# source files for library
LIBsrc_SRC := $(filter-out $(EXEsrc_SRC),$(wildcard $(DIR)/*.cpp))

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
	NM="$(NM)" DLLTOOL="$(DLLTOOL)" \
	CXX="$(FXX)" REALCXX="$(CXX) $(CXXFLAGS)" \
	CC="$(FCC)" REALCC="$(CC) $(CFLAGS)" MATH="$(MATH)" \
	$(MCC) -o $@ $(CPPFLAGS) $(CFLAGS) $^ $(LDLIBS)

ALLDEP += $(LIBsrc_DEP) $(EXEsrc_DEP)
ALLLIB += $(LIBsrc) $(SHAREDLIBsrc)
ALLEXE += $(EXEsrc_EXE)
