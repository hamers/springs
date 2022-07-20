CC=gcc
CXX=g++-11

RM=rm -f
CPPFLAGS=-g -O3
LDFLAGS=-g  
LDLIBS=

SRCS=springs.cpp integrator.cpp tools.cpp parameters.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: springs

springs: $(OBJS)
	$(CXX) $(LDFLAGS) -o springs $(OBJS) $(LDLIBS)

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) springs
