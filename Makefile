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
CPPFLAGS = -gstabs+ -I/u24/local/include $(CFLAGS)
LDFLAGS = -L/u24/local/lib -lgsl -lgslcblas

# PROGRAMS = darwin mrmcpfl c60gradmin plotstrupdf plain2eye eye2plain
PROGRAMS =
TESTS   = $(patsubst %.cpp,%,$(wildcard *Test*.cpp))
SOURCES = BGAlib.cpp BGAlib3d.cpp
HEADERS = $(SOURCES:%.cpp=%.hpp)
OBJECTS = $(SOURCES:%.cpp=%.o)

########################################################################
# most common targets:
all:	 	BGAlib3d.o
# all:	 	$(PROGRAMS)
# all:	 	tests
lib_objects: $(OBJECTS)

########################################################################
# simulations {{{
%.o: %.cpp %.hpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<

#####################################################################}}}
# test programs {{{

tests: $(TESTS)

ssTest01: ssTest01.o $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

molTest01: molTest01.o $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

molTest02: molTest02.o $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

molTest03: molTest03.o $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

molTest04: molTest04.o $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

molTest05: molTest05.o $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

molTest06: molTest06.o $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

molTest07: molTest07.o $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

molTest08: molTest08.o $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

# molTest09: molTest09.o $(OBJECTS) $(HEADERS)
# 	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)
#
# molTest10: molTest10.o $(OBJECTS) $(HEADERS)
# 	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

popTest01: popTest01.o $(OBJECTS) $(HEADERS)
	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

# popTest02: popTest02.o $(OBJECTS) $(HEADERS)
# 	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)
#
# popTest03: popTest03.o $(OBJECTS) $(HEADERS)
# 	$(CC) -o $@ $@.o $(OBJECTS) $(LDFLAGS)

#####################################################################}}}
# utility targets {{{
tags:	$(SOURCES) $(HEADERS) $(PROGRAMS:%=%.cpp)
	ctags $^

clean:
	rm -f $(PROGRAMS) $(OBJECTS) $(PROGRAMS:%=%.o) \
	    a.out \
	    $(TESTS) $(TESTS:%=%.o) \
	    core.[1-9]*

list:
	@printf "PROGRAMS:\n"
	@printf "  %s\n" $(PROGRAMS)
	@printf "HEADERS:\n"
	@printf "  %s\n" $(HEADERS)
	@printf "SOURCES:\n"
	@printf "  %s\n" $(SOURCES)
	@printf "TESTS:\n"
	@printf "  %s\n" $(TESTS:%=%.cpp)

#}}}

# vim: set foldmethod=marker:
