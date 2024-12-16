CC=g++
CFLAGS=-O3
LDFLAGS=#-static

FLORA: FLORA.cpp PDBParser.hpp pstream.h PDBFiller.hpp Superpose.hpp IdealRNA.hpp GeometryTools.hpp cssr.hpp BondLengths.hpp BaseConformation.hpp AtomicClashes.hpp AdjustPosition.hpp MissingRNAatom.hpp SSS.hpp PDBOptimizer.hpp Chirality.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

Arena2: Arena2.cpp PDBParser.hpp pstream.h PDBFiller.hpp Superpose.hpp IdealRNA.hpp GeometryTools.hpp cssr.hpp BondLengths.hpp BaseConformation.hpp AtomicClashes.hpp AdjustPosition.hpp MissingRNAatom.hpp SSS.hpp PDBOptimizer.hpp Chirality.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

clean:
	rm FLORA
