DIR := doc

HTML_OUTPUT_DIR := $(DIR)/html

.PHONY: doc doc-html

doc: doc-html

doc-html: config.h $(LIBsrc_SRC)
	( cat $(DIR)/Doxyfile ; \
	  echo "INPUT = src README" ; \
	  echo "OUTPUT_DIRECTORY = $(HTML_OUTPUT_DIR)" ; \
	  echo "EXCLUDE = $(ALLDEP)"; \
	) | doxygen -

clean::
	-rm -rf $(HTML_OUTPUT_DIR)
