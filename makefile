# compilation flags
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc

#Link Libraries 
LDLIBSOPTIONS=-L sdsl-lite/SDSL_INSTALL_PATH/lib 

#Include headers
INCLUDEOPTIONS=-I sdsl-lite/SDSL_INSTALL_PATH/include

# main executables 
EXECS=circpfp.x circpfpr.x circpfpd.x
# executables not using threads (and therefore not needing the thread library)
EXECS_NT=circpfpNT.x circpfprNT.x circpfpdNT.x bebwtNT.x bebwtNT64.x invertNT.x
#pfebwtNT.x

# targets not producing a file declared phony
.PHONY: all clean install_sdsl

all: $(EXECS) $(EXECS_NT)

gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

csais.o: csais.cpp csais.h
	$(CC) $(CFLAGS) -c -o $@ $< -ldsl -ldivsufsort -ldivsufsort64

csais64.o: csais.cpp csais.h
	$(CC) $(CFLAGS) -c -o $@ $< -ldsl -ldivsufsort -ldivsufsort64 -DP64

#invertNT.o: invertebwt.cpp common.hpp
#	$(CC) $(CFLAGS) -c -o $@ $< -ldsl -ldivsufsort -ldivsufsort64

invertNT.x: invertebwt.cpp common.hpp
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lsdsl -ldivsufsort -ldivsufsort64

circpfp.x: circpfp.cpp circpfp.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ circpfp.cpp malloc_count.o utils.o xerrors.o -ldl -lz -pthread

circpfpNT.x: circpfp.cpp malloc_count.o utils.o 
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lz -ldl -DNOTHREADS

circpfpr.x: circpfpmw.cpp circpfpmw.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ circpfpmw.cpp malloc_count.o utils.o xerrors.o -ldl -lz -pthread

circpfprNT.x: circpfpmw.cpp malloc_count.o utils.o 
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lz -ldl -DNOTHREADS

circpfpd.x: circpfpd.cpp circpfpd.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ circpfpd.cpp malloc_count.o utils.o xerrors.o -ldl -lz -pthread

circpfpdNT.x: circpfpd.cpp malloc_count.o utils.o 
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lz -ldl -DNOTHREADS

bebwtNT.x: ebwt.cpp parse.hpp dictionary.hpp pfp.hpp common.hpp malloc_count.o utils.o gsa/gsacak.o csais.o
	$(CXX) $(CXX_FLAGS) -o $@ ebwt.cpp malloc_count.o gsa/gsacak.o utils.o csais.o -ldl -lsdsl -ldivsufsort -ldivsufsort64

bebwtNT64.x: ebwt.cpp parse.hpp dictionary.hpp pfp.hpp common.hpp malloc_count.o utils.o gsa/gsacak64.o csais64.o
	$(CXX) $(CXX_FLAGS) -o $@ ebwt.cpp malloc_count.o gsa/gsacak64.o utils.o csais64.o -DP64 -DM64 -ldl -lsdsl -ldivsufsort -ldivsufsort64

clean:
	rm -f $(EXECS) $(EXECS_NT) *.o gsa/*.o 