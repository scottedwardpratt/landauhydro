include ../makefile_defs.mk

INCLUDELIBDIRS = -I${LANDAU_INSTALLDIR}/include -I${EIGEN_HOME} -L${LANDAU_INSTALLDIR}/lib

landaulibs : 
	make -C ${LANDAU_HOME} install

iterate : iterate.cc landaulibs
	${CPP} -o iterate iterate.cc ${CFLAGS} ${INCLUDELIBDIRS} -llandau
