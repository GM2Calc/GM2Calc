DIR          := examples
MODNAME      := examples

# source files with main()
EXEexamples_SRC := \
	$(DIR)/example-gm2calc.cpp \
	$(DIR)/example-slha.cpp

EXEexamples_OBJ := $(EXEexamples_SRC:.cpp=.o)

EXEexamples_DEP := $(EXEexamples_OBJ:.o=.d)

EXEexamples_EXE := $(patsubst $(DIR)/%.o, $(BINDIR)/%.x, $(EXEexamples_OBJ))

.PHONY: examples

clean::
	-rm -f $(EXEexamples_DEP)
	-rm -f $(EXEexamples_OBJ)
	-rm -f $(EXEexamples_EXE)

examples: $(EXEexamples_EXE) make.args

$(EXEexamples_DEP) $(EXEexamples_OBJ): \
	override CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS)

$(BINDIR)/%.x: $(DIR)/%.o $(LIBsrc) | $(BINDIR)
	$(CXX) -o $@ $^ $(LDLIBS)

ALLDEP += $(EXEexamples_DEP)
