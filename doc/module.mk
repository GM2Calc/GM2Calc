DIR := doc

HTML_OUTPUT_DIR := $(DIR)/html

.PHONY: doc doc-html

doc: doc-html

doc-html: $(LIBsrc_SRC)
	( cat $(DIR)/Doxyfile ; \
	  echo "PROJECT_NUMBER = $(VERSION)" ; \
	  echo "INPUT = src README $(DIR)/mainpage.dox" ; \
	  echo "OUTPUT_DIRECTORY = $(HTML_OUTPUT_DIR)" ; \
	  echo "EXCLUDE = $(ALLDEP)"; \
	) | doxygen -

clean::
	-rm -rf $(HTML_OUTPUT_DIR)
