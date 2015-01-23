SVNREV=$(shell svn info .. | grep Revision: | cut -d ' ' -f 2 2> /dev/null)
ifeq ($(SVNREV),)
	SVNREV="-1"
endif

all: html

html: Doxyfile
	@doxygen Doxyfile

Doxyfile: Doxyfile.tmpl Doxyfile.defs
	@cat Doxyfile.tmpl Doxyfile.defs > Doxyfile
	@sed -i -e 's/####REVISION####/$(SVNREV)/' Doxyfile

Doxyfile.defs:
	$(error "Error: you have to do-configure BACI before the documentation can be built.")

.PHONY: Doxyfile clean

clean:
	@rm -rf Doxyfile html/ *.tmp
