##################################################################
##################################################################
#		MAKEFILE for WABBIT
##################################################################
##################################################################
all: directories wabbit wabbit-post #python
##################################################################
include LIB/fortran.mk	# includes the makefile of fortran library
include LIB/python.mk 	# includes the makefile of python library
##################################################################

#================================================================
# Information about the makefile
#================================================================
info:
	@echo -e "\n\n\n"
	@echo -e  "command \t \t info"
	@echo -e  "-------------------------------------------------------"
	@echo -e "\e[1m make\e[0m \t\t \t generates all fortran binaries"
	@echo -e "\e[1m make lib\e[0m \t\t fortran wabbit library libwabbit.a"
	@echo -e "\e[1m make wabbit\e[0m \t\t generates wabbit main"
	@echo -e "\e[1m make test\e[0m \t\t performs unit tests"
	@echo -e "\e[1m make doc\e[0m \t \t generates documentation"
	@echo -e "\e[1m make python\e[0m \t \t generate wabbit python lib"
	@echo -e "\e[1m make clean\e[0m \t \t cleans binaries"
	@echo -e "\e[1m make clean-python\e[0m \t clean python binaries"
	@echo -e "\e[1m make clean-fortran\e[0m \t clean fortran binaries"
	@echo -e  "-------------------------------------------------------"
	@echo -e "\n\n\n"


.PHONY: doc test directories
#================================================================
# Library of all wabbit modules
#================================================================
lib: directories libwabbit.a

#================================================================
# Documentation using doxygen
#================================================================
doc:
	doxygen doc/doc_configuration
	firefox doc/output/html/index.html &
#================================================================
# Unit Testing
#================================================================
test: all 
	./TESTING/runtests.py

compression_test: all
	cd ./TESTING/compression; ./compression_test.py --wabbit-dir="../../" -p

#================================================================
# If the object directory doesn't exist, create it.
#================================================================
directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}

#================================================================
# If the object directory doesn't exist, create it.
#================================================================
clean: clean-fortran clean-python

##################################################################
# Remarks:
# 1. the @ forces make to execute the command without printing it
# 	 first
# 2. to execute a make in another directory append the command using
#    semicolon (;)
##################################################################
