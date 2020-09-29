include ../makefile_defs.mk

INCLUDELIBDIRS = -I${LANDAU_INSTALLDIR}/include -I${EIGEN_HOME} -L${LANDAU_INSTALLDIR}/lib

landaulibs : 
	make -C ${LANDAU_HOME} install

hydro : hydromain.cc landaulibs
	${CPP} -o hydro hydromain.cc ${CFLAGS} ${INCLUDELIBDIRS} -llandau


