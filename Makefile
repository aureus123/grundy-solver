
# Makefile

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CC = g++

CCOPT = -O9 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -Wall -std=c++11 -Wno-write-strings
CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio127/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio127/concert
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I.
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
LIBS = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lm -lconcert -lilocplex -lcplex -pthread

all: grundy gengraph genkneser

grundy: grundy.cpp
	$(CC) -o $@ $^ $(CCFLAGS) $(LIBS) $(CCLNFLAGS)

gengraph: gengraph.cpp
	$(CC) -o $@ $^ $(CCFLAGS) -lm

genkneser: genkneser.cpp
	$(CC) -o $@ $^ $(CCFLAGS) -lm

.PHONY: clean

clean:
	rm -f grundy gengraph genkneser

