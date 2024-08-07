CXXFLAGS=-DNDEBUG -std=c++11 -O3 -w -fPIC
CC=g++
#CC=x86_64-w64-mingw32-g++
LDFLAGS=-lstdc++ #-static

CSSR: cssr_struct.o PDBParser.o CSSR.o
	${CC} -o CSSR ${CXXFLAGS} cssr_struct.o PDBParser.o CSSR.o

CSSR.o: CSSR.cpp
	${CC} -c -o $@ ${CXXFLAGS} CSSR.cpp
PDBParser.o: PDBParser.cpp
	${CC} -c -o $@ ${CXXFLAGS} PDBParser.cpp
cssr_struct.o: cssr_struct.cpp
	${CC} -c -o $@ ${CXXFLAGS} cssr_struct.cpp

clean:
	rm *.o CSSR
