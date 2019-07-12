#######################################################################################
#                                                                                     #
# Author: Paolo Paronuzzi        Unibo     09-07-2018	                              #
#                                                                                     #
#######################################################################################
#    You can modify this example makefile to fit for your own program. 	 	      #
#    Usually, you only need to change the five CHANGEME entries below. 		      #
#######################################################################################

# CHANGEME: This should be the name of your executable
EXE = KSEP

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS = main.o global_variables.o Compact_OK.o Graph_v4.o global_functions.o StdPath_TerminalModel.o StdPath_withLongPath.o BilevelModel.o BilevelModel_withLeaf.o SmartModel.o LP_Model.o CLIQUE.o

# Linux and CPLEX VERSION (LIBS and INCLUDE files)
SYSTEM = x86-64_linux
LIBFORMAT = static_pic
CPLEXDIR = /opt/ibm/ILOG/CPLEX_Studio1271

CPLEXLIBDIR = $(CPLEXDIR)/cplex/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CPLEXDIR)/concert/lib/$(SYSTEM)/$(LIBFORMAT)
LP_CPLEX_INCLUDE= $(CPLEXDIR)/cplex/include
LP_CONCERT_INCLUDE = $(CPLEXDIR)/concert/include
INCDIR = -I. -I$(LP_CONCERT_INCLUDE) -I$(LP_CPLEX_INCLUDE)
	
# CHANGEME: Additional libraries
ADDLIBS = -L$(CONCERTLIBDIR) -lconcert -L$(CPLEXLIBDIR) -lilocplex -lcplex -lm -lpthread 

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = $(INCDIR)

# CHANGEME: Directory to the sources for the (example) problem definition
# files
SRCDIR = src
VPATH = src
OBJDIR = obj

#######################################################################################
#  Usually, you don't have to change anything below.                                  #
#######################################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O3 -pipe -DIL_STD -DNDEBUG -pedantic-errors -DBONMIN_BUILD -Wno-cast-qual -Wno-long-long -Wno-unknown-pragmas -Wno-unused-result

INCL += $(ADDINCFLAGS)

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $(OBJDIR)/$$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) $(INCL) -o $@ $$bla $(LIBS) $(ADDLIBS)	
	
clean:
	rm -rf $(EXE) $(OBJDIR)/$(OBJS)

remake:
	make clean;
	make;
	
.cpp.o:
	@mkdir -p $(OBJDIR);\
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $(OBJDIR)/$@ `test -f '$<' || echo '$(SRCDIR)/'`$<

.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $(OBJDIR)/$@ `if test -f '$<'; \
	then $(CYGPATH_W) '$<';	else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`
