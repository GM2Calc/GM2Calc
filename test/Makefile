DIR        := .
SRC_DIR    := $(DIR)/../src

CPPFLAGS   := -I$(SRC_DIR)
CXXFLAGS   := -O2 -std=c++11
BOOSTFLAGS := -I/usr/include
EIGENFLAGS := -I/usr/include/eigen3
MATH       := math

TEST_M     := $(DIR)/test_mathlink_interface.m
TEST_SRC   := $(DIR)/test_c_interface.cpp
TEST_SH    := $(DIR)/test_c_examples.sh $(DIR)/test_points.sh
TEST_DEP   := $(TEST_SRC:.cpp=.d)
TEST_OBJ   := $(TEST_SRC:.cpp=.o)
TEST_EXE   := $(TEST_SRC:.cpp=.x)
TEST_LOG   := \
	$(patsubst %.m, %.m.log, $(filter %.m, $(TEST_M))) \
	$(patsubst %.cpp, %.x.log, $(filter %.cpp, $(TEST_SRC))) \
	$(patsubst %.sh, %.sh.log, $(filter %.sh, $(TEST_SH)))

LIBGM2Calc := $(SRC_DIR)/libgm2calc.a

all: $(TEST_EXE)

ifneq ($(MAKECMDGOALS),clean)
-include $(TEST_DEP)
endif

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

%.m.log: %.m $(TEST_M)
	@rm -f $@
	@echo "**************************************************" >> $@;
	@echo "* executing test: $< " >> $@;
	@echo "**************************************************" >> $@;
	@"$(MATH)" -run "Get[\"$<\"];" >> $@ 2>&1; \
	if [ $$? = 0 ]; then echo "$<: OK"; else echo "$<: FAILED"; fi

%.log: %
	@rm -f $@
	@echo "**************************************************" >> $@;
	@echo "* executing test: $< " >> $@;
	@echo "**************************************************" >> $@;
	@./$< > $@ 2>&1; \
	if [ $$? = 0 ]; then echo "$<: OK"; else echo "$<: FAILED"; fi
