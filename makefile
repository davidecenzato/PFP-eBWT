# compilation flags
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc

#Link Libraries 
LDLIBSOPTIONS=-L sdsl-lite/SDSL_INSTALL_PATH/lib 

#Include headers
INCLUDEOPTIONS=-I sdsl-lite/SDSL_INSTALL_PATH/include

# main executables 
EXECS=newscan.x 
# executables not using threads (and therefore not needing the thread library)
EXECS_NT=newscanNT.x pfebwtNT.x 
#pfebwtNT.x

# targets not producing a file declared phony
.PHONY: all clean install_sdsl

all: $(EXECS) $(EXECS_NT)

gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

newscan.x: newscan.cpp newscan.hpp malloc_count.o utils.o xerrors.o 
	$(CXX) $(CXX_FLAGS) -o $@ newscan.cpp malloc_count.o utils.o xerrors.o -ldl -lz -pthread

newscanNT.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -lz -ldl -DNOTHREADS

pfebwtNT.x: ebwt.cpp parse.hpp dictionary.hpp pfp.hpp common.hpp malloc_count.o utils.o gsa/gsacak.o
	$(CXX) $(CXX_FLAGS) -o $@ ebwt.cpp malloc_count.o gsa/gsacak.o utils.o ${INCLUDEOPTIONS} ${LDLIBSOPTIONS} -ldl -lsdsl -ldivsufsort -ldivsufsort64

clean:
	rm -f $(EXECS) $(EXECS_NT) *.o gsa/*.o

install_sdsl:
	./sdsl-lite/install.sh

unistall_sdsl:
	./sdsl-lite/uninstall.sh