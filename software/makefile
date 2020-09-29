include ../makefile_defs.mk

install : installdirs ${LANDAU_INSTALLDIR}/lib/liblandau.a ${LANDAU_INSTALLDIR}/include/defs.h ${LANDAU_INSTALLDIR}/include/randy.h ${LANDAU_INSTALLDIR}/include/landau.h ${LANDAU_INSTALLDIR}/include/constants.h ${LANDAU_INSTALLDIR}/include/parametermap.h ${LANDAU_INSTALLDIR}/include/eos.h

installdirs : 
	mkdir -p ${LANDAU_INSTALLDIR}/include;\
	mkdir -p ${LANDAU_INSTALLDIR}/lib;\

clean :
	rm -f build/*.o lib/*.a

INC=-I include -I ${EIGEN_HOME}

uninstall :
	rm -f ${LANDAU_INSTALLDIR}/lib/liblandau.a ${LANDAU_INSTALLDIR}/include/defs.h ${LANDAU_INSTALLDIR}/include/parametermap.h ${LANDAU_INSTALLDIR}/include/randy.h ${LANDAU_INSTALLDIR}/include/landau.h ${LANDAU_INSTALLDIR}/include/eos.h 

LANDAU_HFILES = include/landau.h include/defs.h include/parametermap.h include/randy.h include/constants.h include/eos.h

LANDAU_OBJFILES = build/landau.o build/cell.o build/mesh.o build/randy.o build/parametermap.o build/eos.o

####################

${LANDAU_INSTALLDIR}/include/landau.h : include/landau.h
	cp -f include/landau.h ${LANDAU_INSTALLDIR}/include/
	
${LANDAU_INSTALLDIR}/include/defs.h : include/defs.h
	cp -f include/defs.h ${LANDAU_INSTALLDIR}/include/

${LANDAU_INSTALLDIR}/include/constants.h : include/constants.h
	cp -f include/constants.h ${LANDAU_INSTALLDIR}/include/

${LANDAU_INSTALLDIR}/include/randy.h : include/randy.h
	cp -f include/randy.h ${LANDAU_INSTALLDIR}/include/

${LANDAU_INSTALLDIR}/include/parametermap.h : include/parametermap.h
	cp -f include/parametermap.h ${LANDAU_INSTALLDIR}/include/

${LANDAU_INSTALLDIR}/include/eos.h : include/eos.h
	cp -f include/eos.h ${LANDAU_INSTALLDIR}/include/

${LANDAU_INSTALLDIR}/lib/liblandau.a : lib/liblandau.a
	cp -r lib/liblandau.a ${LANDAU_INSTALLDIR}/lib/

######################

lib/liblandau.a : ${LANDAU_HFILES} ${LANDAU_OBJFILES}
	rm -f lib/liblandau.a;\
	ar -ru lib/liblandau.a ${LANDAU_OBJFILES};
	
build/landau.o : src/landau.cc ${LANDAU_HFILES}
	${CPP} -c ${CFLAGS} ${INC} -o build/landau.o src/landau.cc

build/cell.o : src/cell.cc ${LANDAU_HFILES}
	${CPP} -c ${CFLAGS} ${INC} -o build/cell.o src/cell.cc

build/mesh.o : src/mesh.cc ${LANDAU_HFILES}
	${CPP} -c ${CFLAGS} ${INC} -o build/mesh.o src/mesh.cc

build/randy.o : src/randy.cc include/randy.h
	${CPP} -c ${CFLAGS} ${INC} -o build/randy.o src/randy.cc

build/parametermap.o : src/parametermap.cc include/parametermap.h
	${CPP} -c ${CFLAGS} ${INC} -o build/parametermap.o src/parametermap.cc

build/eos.o : src/parametermap.cc src/eos.cc include/eos.h
	${CPP} -c ${CFLAGS} ${INC} -o build/eos.o src/eos.cc
	
############

