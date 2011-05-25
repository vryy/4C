
baci:
	cd ../.. && $(MAKE) `cat executable.txt 2> /dev/null || grep ^baci-[a-z]*: Makefile | awk '{ print substr($$1,0,length($$1)-1) }'`
