###################################################################
# This Makefile was created using the bat-project script
# for project binomial
# bat-project is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://mpp.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CXXFLAGS and LIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# List of all class (model) sources used in the program,
# separated by spaces. A backslash indicates continuation
# on the next line
CXXSRCS = RefHandler.cxx FitHandler.cxx LinAlg.cxx ConfigHandler.cxx WeightsHandler.cxx Counter.cxx Power.cxx OutputHandler.cxx
#datareading.cxx
# List of all program sources used in the program,
# separated by spaces. A backslash indicates continuation
# on the next line
PRGSRCS = qr.cxx 

# compiler and flags
CXX       = g++
CXXFLAGS  = -I /home/iwsatlas1/fabiankl/pythiatune/local/Eigen -I /home/iwsatlas1/fabiankl/pythiatune/build/root-6.06.04/include -pg -g -O2 -Wall -fPIC -Wno-deprecated 
LD        = /usr/bin/ld -m elf_x86_64
LDFLAGS   = -g -O2 

# ----------------------------------------------------------------------
# The following definitions rely on the script bat-config being
# available in $PATH. If BAT is not installed in the standard system
# directories, update $PATH accordingly.

CXXFLAGS += -fopenmp -std=c++11 -Wno-deprecated-declarations -m64
LIBS := -L/home/iwsatlas1/fabiankl/pythiatune/build/root-6.06.04/lib -fopenmp -lCore -lMatrix


# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS = $(addsuffix .o,$(basename $(CXXSRCS)))
MYPROGS = $(basename $(PRGSRCS))
PRGOBJS = $(addsuffix .o,$(basename $(PRGSRCS)))

GARBAGE = $(CXXOBJS) $(PRGOBJS) link.d $(MYPROGS)

# targets
all : $(MYPROGS)

.PHONY : all clean print

link.d : $(addsuffix .h,$(basename $(CXXSRCS))) $(CXXSRCS) $(PRGSRCS)
	$(CXX) -MM $(CXXFLAGS) $(filter-out %.h,$^) > link.d;
	@$(foreach prog,$(MYPROGS), echo $(prog) : $(prog).o >> link.d;)

-include link.d

$(CXXOBJS) $(PRGOBJS) :
	$(CXX) $(CXXFLAGS) -c $(filter $(basename $@).%,$(filter-out %.h,$^)) -o $@

$(MYPROGS) : $(CXXOBJS)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

clean :
	rm -f $(GARBAGE)

print :
	@echo compiler  : $(CXX)
	@echo c++ srcs  : $(CXXSRCS) $(PRGSRCS)
	@echo c++ objs  : $(CXXOBJS) $(PRGOBJS)
	@echo c++ flags : $(CXXFLAGS)
	@echo ld flags  : $(LDFLAGS)
	@echo libs      : $(LIBS)
