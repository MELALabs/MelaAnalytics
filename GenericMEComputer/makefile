COMPILEPATH          = $(PWD)/
BASEINCLUDE          = $(COMPILEPATH)/../../
MELADIR              = $(BASEINCLUDE)ZZMatrixElement/MELA
INCLUDEDIR           = $(COMPILEPATH)interface/
SRCDIR               = $(COMPILEPATH)src/
BINDIR               = $(COMPILEPATH)bin/
LIBDIR               = $(COMPILEPATH)lib/
TESTDIR              = $(COMPILEPATH)test/
RUNDIR               = $(COMPILEPATH)
LIB                  = libMelaAnalyticsGenericMEComputer.so

ROOTCFLAGS    = $(shell root-config --cflags) -Lrootlib -lEG
ROOTGLIBS     = $(shell root-config --glibs) -Lrootlib -lEG

ARCH         := $(shell root-config --arch)

CXX           = g++
CXXFLAGS      = -fPIC -g -O2 $(ROOTCFLAGS) -I$(ROOFITSYS)/include/ -I$(INCLUDEDIR) -I$(MELADIR)/interface/
LINKERFLAGS   = -Wl,-rpath=$(LIBDIR),-soname,$(LIB)

NGLIBS         = $(ROOTGLIBS)
# Hack here, because RooFit is removed from ROOT:
NGLIBS        += -L $(ROOFITSYS)/lib/ -lMinuit -lRooFitCore -lRooFit #-lboost_regex
# Libraries for common user packages
NGLIBS        += -lZZMatrixElementMELA
# Filter out libNew because it leads to floating-point exceptions in versions of ROOT prior to 6.08.02
# See the thread https://root-forum.cern.ch/t/linking-new-library-leads-to-a-floating-point-exception-at-startup/22404
GLIBS          = $(filter-out -lNew, $(NGLIBS))


SOURCES = $(wildcard $(SRCDIR)*.cc)
OBJECTS = $(SOURCES:.cc=.o)

.PHONY: all help compile clean

all: $(LIB)


# General rule for making object files
%.d:	$(SRCDIR)%.cc $(BINDIR)%.cc $(SRCDIR)%.cxx $(BINDIR)%.cxx
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CXX) -MM -MT $@ -MT ${@:.d=.o} $(CXXFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	$(SRCDIR)%.cc $(SRCDIR)%.cxx $(BINDIR)%.cc $(BINDIR)%.cxx
	$(QUIET) echo "Compiling $<"; \
	$(CXX) $(CXXFLAGS) $< -c -o $@ $(GLIBS)


$(LIB):	$(OBJECTS)
	$(QUIET) echo "Linking $(LIB)"; \
	$(CXX) $(LINKERFLAGS) -shared $(OBJECTS) -o $@; \
        mkdir -p $(LIBDIR); \
	mv $(LIB) $(LIBDIR)/$(LIB)

clean:
	rm -rf $(LIBDIR)
	rm -f $(SRCDIR)*.o
	rm -f $(SRCDIR)*.so
	rm -f $(SRCDIR)*.d
	rm -f $(BINDIR)*.o
	rm -f $(BINDIR)*.so
	rm -f $(BINDIR)*.d
	rm -rf $(RUNDIR)Pdfdata
	rm -f $(RUNDIR)input.DAT
	rm -f $(RUNDIR)process.DAT
	rm -f $(RUNDIR)ffwarn.dat
	rm -f $(RUNDIR)br.sm*
	rm -f $(RUNDIR)*.cc
	rm -f $(RUNDIR)*.o
	rm -f $(RUNDIR)*.so
	rm -f $(RUNDIR)*.d
	rm -f $(RUNDIR)*.pcm
	rm -f $(RUNDIR)*.py
	rm -f $(RUNDIR)*.pyc
	rm -rf $(TESTDIR)Pdfdata
	rm -f $(TESTDIR)input.DAT
	rm -f $(TESTDIR)process.DAT
	rm -f $(TESTDIR)ffwarn.dat
	rm -f $(TESTDIR)br.sm*
	rm -f $(TESTDIR)*.o
	rm -f $(TESTDIR)*.so
	rm -f $(TESTDIR)*.d
	rm -f $(TESTDIR)*.pcm
	rm -f $(TESTDIR)*.py
	rm -f $(TESTDIR)*.pyc
