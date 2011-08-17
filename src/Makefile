# Name of the executable
EXE=minimac
OMP_EXE=minimac-omp
########################
# The Files:
TOOLBASE = HaplotypeClipper HaplotypeSet ImputationStatistics MarkovModel MarkovParameters
SRCONLY = Main.cpp
HDRONLY = 

########################
# Handle openmp
ifeq ($(MAKECMDGOALS),openmp)
  OBJDIR_OPT = $(OBJDIR)/omp
  EXE = $(OMP_EXE)
endif

########################
# Include the base Makefile
PARENT_MAKE = Makefile.src
include ../Makefile.inc

########################
# Handle openmp
openmp: USER_COMPILE_VARS = -fopenmp
openmp: opt
	echo $(EXE)

# Also need to build openmp
all:
	$(MAKE) --no-print-directory openmp
########################
# Handle openmp
USER_REMOVES += -rm -rf $(OBJDIR)/omp/*.o $(BINDIR)/$(OMP_EXE)
