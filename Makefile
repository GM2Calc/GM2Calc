# Package information
PKGNAME         := gm2calc
VERSION         := 0.2.10

# Variables for compilation
CXX             := g++
CPPFLAGS        := -Isrc
CXXFLAGS        := -O2 -std=c++11
CXX_DEP_GEN     := g++
FC              := gfortran
FFLAGS          := -O2 -frecursive
FLIBS           := -lgfortran -lm
FOR_DEP_GEN     := gfortran
MAKELIB         := ar cru
BLASLIBS        := -lblas
BOOSTFLAGS      := -I/usr/include
EIGENFLAGS      := -I/usr/include/eigen3
LAPACKLIBS      := -llapack
LIBEXT          := .a
BINDIR          := bin
CONFIG_H        := src/config.h

.PHONY:         all allexec alllib clean depend make.args

all: alllib allexec make.args

clean::
	-rm -f $(CONFIG_H)
	-rm -rf $(BINDIR)
	-rm -f make.args

include src/module.mk
include doc/module.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),depend)
ifneq ($(MAKECMDGOALS),doc)
ifneq ($(MAKECMDGOALS),make.args)
ifeq ($(findstring doc-,$(MAKECMDGOALS)),)
ifeq ($(findstring print-,$(MAKECMDGOALS)),)
-include $(ALLDEP)
endif
endif
endif
endif
endif
endif

allexec:  $(ALLEXE)
alllib:   $(ALLLIB)
depend:   $(ALLDEP)

$(BINDIR):
	mkdir $(BINDIR)

$(CONFIG_H): Makefile
	rm -f $@-t $@
	{ echo '/* DO NOT EDIT! GENERATED AUTOMATICALLY! */'; \
	  echo ''; \
	  echo '/* program version */'; \
	  echo '#define GM2CALC_VERSION "$(VERSION)"'; \
	} | sed '/""/d' > $@-t
	mv $@-t $@

make.args:
	rm -f $@-t $@
	{ echo 'CXX="$(CXX)"' \
	       'CPPFLAGS="$(CPPFLAGS)"' \
	       'CXXFLAGS="$(CXXFLAGS)"' \
	       'CXX_DEP_GEN="$(CXX_DEP_GEN)"' \
	       'FC="$(FC)"' \
	       'FFLAGS="$(FFLAGS)"' \
	       'FLIBS="$(FLIBS)"' \
	       'FOR_DEP_GEN="$(FOR_DEP_GEN)"' \
	       'MAKELIB="$(MAKELIB)"' \
	       'BLASLIBS="$(BLASLIBS)"' \
	       'BOOSTFLAGS="$(BOOSTFLAGS)"' \
	       'EIGENFLAGS="$(EIGENFLAGS)"' \
	       'LAPACKLIBS="$(LAPACKLIBS)"'; \
	} > $@-t
	mv $@-t $@

%.d: %.cpp | $(CONFIG_H)
	$(CXX_DEP_GEN) $(CPPFLAGS) -MM -MP -MG -o $@ -MT '$*.o' $^

%.d: %.f | $(CONFIG_H)
	$(FOR_DEP_GEN) $(CPPFLAGS) -cpp -MM -MP -MG $^ -MT '$*.o' | \
	sed 's|.*\.o:|$*.o:|' > $@

print-% : ; @echo $* = $($*)

tag:
	git tag v$(VERSION)

release-tag:
	git archive --worktree-attributes --prefix=gm2calc-v$(VERSION)/ \
		--output=gm2calc-v$(VERSION).tar.gz v$(VERSION)
	md5sum gm2calc-v$(VERSION).tar.gz > gm2calc-v$(VERSION).tar.gz.md5

release-head:
	$(eval GIT_HEAD_DESCR := $(shell git describe --tags HEAD))
	git archive --worktree-attributes --prefix=gm2calc-$(GIT_HEAD_DESCR)/ \
		--output=gm2calc-$(GIT_HEAD_DESCR).tar.gz v$(GIT_HEAD_DESCR)
	md5sum gm2calc-$(GIT_HEAD_DESCR).tar.gz \
		> gm2calc-$(GIT_HEAD_DESCR).tar.gz.md5
