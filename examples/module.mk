DIR          := examples
MODNAME      := examples

# source files with main()
EXEexamples_SRC := \
	$(DIR)/example-gm2calc_c.c \
	$(DIR)/example-gm2calc.cpp \
	$(DIR)/example-slha_c.cpp \
	$(DIR)/example-slha.cpp

EXEexamples_OBJ := \
	$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEexamples_SRC))) \
	$(patsubst %.c, %.o, $(filter %.c, $(EXEexamples_SRC)))

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

$(BINDIR)/%_c.x: $(DIR)/%_c.o $(LIBsrc) | $(BINDIR)
	$(CC) -o $@ $^ $(LDLIBS) $(CLIBS)

$(BINDIR)/%.x: $(DIR)/%.o $(LIBsrc) | $(BINDIR)
	$(CXX) -o $@ $^ $(LDLIBS)

ALLDEP += $(EXEexamples_DEP)
