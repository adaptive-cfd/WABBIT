##################################################################
##################################################################
#		MAKEFILE for WABBIT
##################################################################
##################################################################
include LIB/fortran.mk		# includes the makefile of fortran
include LIB/python.mk 		# includes the makefile of python
##################################################################

#================================================================
# Information about the makefile
#================================================================
info:
	@echo -e "\n\n\n"
	@echo -e  "command \t \t info"
	@echo -e  "-------------------------------------------------------"
	@echo -e "\e[1m make\e[0m \t\t \t generates all fortran binaries"
	@echo -e "\e[1m make test\e[0m \t\t performs unit tests"
	@echo -e "\e[1m make clean\e[0m \t \t cleans fortran binaries"
	@echo -e "\e[1m make doc\e[0m \t \t generates documentation"
	@echo -e "\e[1m make python\e[0m \t \t generate wabbit python lib"
	@echo -e "\e[1m make python-clean\e[0m \t clean python binaries"
	@echo -e "\e[1m make ctags\e[0m \t generates ctags for your editor"
	@echo -e  "-------------------------------------------------------"
	@echo -e "\n\n\n"

.PHONY: doc test
#================================================================
# Documentation using doxygen
#================================================================
doc:
	doxygen doc/doc_configuration
	firefox doc/output/html/index.html &
#================================================================
# Unit Testing
#================================================================
test:
	./TESTING/runtests.sh

#================================================================
# Generate Ctags from WABBIT source Library
#================================================================
ctags:
	mkdir -p OBJ
	ctags -R --fortran-kinds=+i -f OBJ/wabbit_tags LIB/*
	sed -i -e 's/\.\.\///' OBJ/wabbit_tags


##################################################################
# Remarks:
# 1. the @ forces make to execute the command without printing it
# 	 first
# 2. to execute a make in another directory append the command using
#    semicolon (;)
##################################################################
