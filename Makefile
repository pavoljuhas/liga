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
%.o: %.cpp %.h
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<

#####################################################################}}}
# test programs {{{

test: $(TESTS)

ssTest1: ssTest1.o $(OBJECTS)
	$(CC) -o $@ $^ $(LDFLAGS) 

molTest1: molTest1.o $(OBJECTS)
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
########################################################################
# Here is what people have been up to: {{{
#
# $Log$
# Revision 1.4  2005/01/25 17:35:48  juhas
# *** empty log message ***
#
# Revision 1.3  2005/01/25 17:32:56  juhas
# added test program targets
#
# Revision 1.2  2005/01/25 16:18:31  juhas
# BGAlib.h renamed to BGAlib.hpp
#
# Revision 1.1  2005/01/24 20:24:13  juhas
# *** empty log message ***
#
# vim: set foldmethod=marker: ########################################}}}
