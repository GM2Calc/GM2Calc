DIR        := .
SRC_DIR    := $(DIR)/../src

CPPFLAGS   := -I$(SRC_DIR)
CXXFLAGS   := -O2 -std=c++11
BOOSTFLAGS := -I/usr/include
EIGENFLAGS := -I/usr/include/eigen3

TEST_SRC   := $(DIR)/test_c_interface.cpp
TEST_DEP   := $(TEST_SRC:.cpp=.d)
TEST_OBJ   := $(TEST_SRC:.cpp=.o)
TEST_EXE   := $(TEST_SRC:.cpp=.x)
TEST_LOG   := $(TEST_SRC:.cpp=.x.log)

LIBGM2Calc := $(SRC_DIR)/libgm2calc.a

all: $(TEST_EXE)

-include $(TEST_DEP)

.PHONY: all clean execute-tests

clean:
	-rm -f $(TEST_DEP)
	-rm -f $(TEST_OBJ)
	-rm -f $(TEST_EXE)
	-rm -f $(TEST_LOG)

execute-tests: $(TEST_LOG)

print-% : ; @echo $* = $($*)

$(TEST_DEP) $(TEST_OBJ): \
	override CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS)

%.d: %.cpp
	$(CXX) $(CPPFLAGS) -MM -MP -MG -o $@ -MT '$*.o' $^

%.x: %.o $(LIBGM2Calc)
	$(CXX) -o $@ $^ $(LIBGM2Calc) $(LDLIBS)

%.x.log: %.x
	@rm -f $@
	@echo "**************************************************" >> $@;
	@echo "* executing test: $< " >> $@;
	@echo "**************************************************" >> $@;
	@./$< > $@ 2>&1; \
	if [ $$? = 0 ]; then echo "$<: OK"; else echo "$<: FAILED"; fi
