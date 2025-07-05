# Makefile

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CC = g++ -std=c++20

CCOPT = -O3 -fPIC -fno-strict-aliasing -fexceptions -Wall -g
CCINFLAGS =  -I include/ -I hglib/include/ -I hglib/lib/filtered_vector/include/  -I exactcolors/

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio2211/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio2211/concert
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CPLEXFLAGS = $(CPLEXINFLAGS) $(CPLEXLIBFLAGS) $(CPLEXLNFLAGS)
CPLEXINFLAGS = -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)
CPLEXLIBFLAGS = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)  
CPLEXLNFLAGS = -lconcert -lilocplex -lcplex -ldl 


all: run_tests

stats.o: src/stats.cpp include/stats.hpp
	$(CC) -c -o $@ $< $(CCOPT) $(CCINFLAGS)

graph.o: src/graph.cpp include/graph.hpp include/params.hpp
	$(CC) -c -o $@ $< $(CCOPT) $(CCINFLAGS)

heur.o: src/heur.cpp include/heur.hpp include/col.hpp include/graph.hpp include/random.hpp
	$(CC) -c -o $@ $< $(CCOPT) $(CCINFLAGS)

compact_ilp.o: src/compact_ilp.cpp include/compact_ilp.hpp \
include/graph.hpp include/stats.hpp include/col.hpp
	$(CC) -c -o $@ $< $(CCOPT) $(CCINFLAGS) $(CPLEXFLAGS)

pricing.o: src/pricing.cpp include/pricing.hpp include/graph.hpp include/stats.hpp include/random.hpp
	$(CC) -c -o $@ $< $(CCOPT) $(CCINFLAGS) $(CPLEXFLAGS)

col.o: src/col.cpp include/col.hpp include/graph.hpp
	$(CC) -c -o $@ $< $(CCOPT) $(CCINFLAGS)

lp.o: src/lp.cpp include/lp.hpp include/col.hpp include/cplex_env.hpp \
include/graph.hpp include/stats.hpp include/params.hpp
	$(CC) -c -o $@ $< $(CCOPT) $(CCINFLAGS) $(CPLEXFLAGS)

main.o: main.cpp include/bp.hpp include/lp.hpp include/graph.hpp \
include/col.hpp include/compact_ilp.hpp include/stats.hpp
	$(CC) -c -o $@ $< $(CCOPT) $(CCINFLAGS) $(CPLEXFLAGS)

tests.o: tests.cpp include/bp.hpp include/lp.hpp include/graph.hpp \
include/col.hpp include/compact_ilp.hpp include/stats.hpp include/params.hpp \
include/pricing.hpp include/heur.hpp
	$(CC) -c -o $@ $< $(CCOPT) $(CCINFLAGS) $(CPLEXFLAGS)

tests_lp.o: tests.cpp include/bp.hpp include/lp.hpp include/graph.hpp \
include/col.hpp include/compact_ilp.hpp include/stats.hpp \
include/pricing.hpp include/heur.hpp
	$(CC) -c -o $@ $< $(CCOPT) $(CCINFLAGS) $(CPLEXFLAGS) -DONLY_RELAXATION

run_tests: tests.o lp.o col.o compact_ilp.o graph.o stats.o pricing.o heur.o\
exactcolors/color.o exactcolors/color_version.h exactcolors/util.o \
exactcolors/rounding_mode.o exactcolors/cliq_enum.o exactcolors/color_parms.o \
exactcolors/graph.o exactcolors/lpcplex.o exactcolors/bbsafe.o \
exactcolors/mwis_grdy.o exactcolors/heap.o \
exactcolors/mwis.o exactcolors/mwis_sewell/mwss_ext.o \
exactcolors/mwis_sewell/wstable.o exactcolors/color_backup.o exactcolors/greedy.o
	$(CC) -o $@ $^ $(CCOPT) $(CCINFLAGS) $(CPLEXFLAGS)

run_tests_lp: tests_lp.o lp.o col.o compact_ilp.o graph.o stats.o pricing.o heur.o\
exactcolors/color.o exactcolors/color_version.h exactcolors/util.o \
exactcolors/rounding_mode.o exactcolors/cliq_enum.o exactcolors/color_parms.o \
exactcolors/graph.o exactcolors/lpcplex.o exactcolors/bbsafe.o \
exactcolors/mwis_grdy.o exactcolors/heap.o \
exactcolors/mwis.o exactcolors/mwis_sewell/mwss_ext.o \
exactcolors/mwis_sewell/wstable.o exactcolors/color_backup.o exactcolors/greedy.o
	$(CC) -o $@ $^ $(CCOPT) $(CCINFLAGS) $(CPLEXFLAGS)

bp: main.o lp.o col.o compact_ilp.o graph.o stats.o \
exactcolors/color.o exactcolors/color_version.h exactcolors/util.o \
exactcolors/rounding_mode.o exactcolors/cliq_enum.o exactcolors/color_parms.o \
exactcolors/graph.o exactcolors/lpcplex.o exactcolors/bbsafe.o \
exactcolors/mwis_grdy.o exactcolors/heap.o \
exactcolors/mwis.o exactcolors/mwis_sewell/mwss_ext.o \
exactcolors/mwis_sewell/wstable.o exactcolors/color_backup.o exactcolors/greedy.o
	$(CC) -o $@ $^ $(CCOPT) $(CCINFLAGS) $(CPLEXFLAGS)