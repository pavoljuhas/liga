########################################################################
# Short Title: make rules for Biosphere Genetic Algorithm
#
# Comments:
#
# $Id$
########################################################################

########################################################################
# compiler flags and file lists:
CC = g++
CPPFLAGS = -g -I/u24/local/include $(CFLAGS)
LDFLAGS = -L/u24/local/lib -lgsl -lgslcblas

# PROGRAMS = darwin mrmcpfl c60gradmin plotstrupdf plain2eye eye2plain
PROGRAMS =
TESTS   = $(patsubst %.cpp,%,$(wildcard *Test*.cpp))
SOURCES = BGAlib.cpp
HEADERS = $(SOURCES:%.cpp=%.hpp)
OBJECTS = $(SOURCES:%.cpp=%.o)

########################################################################
# most common targets:
all:	 	BGAlib.o
# all:	 	$(PROGRAMS)
# all:	 	test
lib_objects: $(OBJECTS)

########################################################################
# simulations {{{
%.o: %.cpp %.hpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<

#####################################################################}}}
# test programs {{{

test: $(TESTS)

ssTest01: ssTest01.o $(OBJECTS)
	$(CC) -o $@ $^ $(LDFLAGS) 

molTest01: molTest01.o $(OBJECTS)
	$(CC) -o $@ $^ $(LDFLAGS) 

molTest02: molTest02.o $(OBJECTS)
	$(CC) -o $@ $^ $(LDFLAGS) 

#####################################################################}}}
# utility targets {{{
tags:	$(SOURCES) $(HEADERS) $(PROGRAMS:%=%.cpp)
	ctags $^

clean:
	rm -f $(PROGRAMS) $(OBJECTS) $(PROGRAMS:%=%.o) \
	    a.out \
	    ssTest*.o \
	    molTest*.o


list:
	@printf "%s\n" $(PROGRAMS:%=%.cpp) $(HEADERS) $(SOURCES)

#}}}

# vim: set foldmethod=marker: 
