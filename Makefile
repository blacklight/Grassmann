SRCDIR=src
INCLUDEDIR=include
PREFIX=/usr
LIB=grassmann
CC=g++
CFLAGS=-Wall -pedantic -pedantic-errors

all:
	${CC} -I${INCLUDEDIR} ${CFLAGS} -fPIC -g -c ${SRCDIR}/matrix.cpp
	${CC} -I${INCLUDEDIR} ${CFLAGS} -fPIC -g -c ${SRCDIR}/vector.cpp
	${CC} -I${INCLUDEDIR} ${CFLAGS} -fPIC -g -c ${SRCDIR}/interpolate.cpp
	${CC} -shared -Wl,-soname,lib$(LIB).so.0 -o lib${LIB}.so.0.0.0 vector.o matrix.o interpolate.o
	ar rcs lib${LIB}.a vector.o matrix.o interpolate.o

install:
	mkdir -p ${PREFIX}/lib
	mkdir -p ${PREFIX}/${INCLUDEDIR}
	mkdir -p ${PREFIX}/share/${LIB}
	install -m 0755 lib${LIB}.so.0.0.0 ${PREFIX}/lib/lib${LIB}.so.0.0.0
	install -m 0644 lib${LIB}.a ${PREFIX}/lib/lib${LIB}.a
	install -m 0644 ${INCLUDEDIR}/${LIB}.hpp ${PREFIX}/${INCLUDEDIR}
	install -m 0644 ${INCLUDEDIR}/${LIB}_exception.hpp ${PREFIX}/${INCLUDEDIR}
	install -m 0644 README ${PREFIX}/share/${LIB}
	install -m 0644 INSTALL ${PREFIX}/share/${LIB}
	install -m 0644 LICENCE ${PREFIX}/share/${LIB}
	install -m 0644 ChangeLog ${PREFIX}/share/${LIB}
	cp -r doc ${PREFIX}/share/${LIB}
	cp -r examples ${PREFIX}/share/${LIB}
	ln -sf ${PREFIX}/lib/lib${LIB}.so.0.0.0 ${PREFIX}/lib/lib${LIB}.so.0

clean:
	rm *.o
	rm lib${LIB}.so.0.0.0
	rm lib${LIB}.a

uninstall:
	rm ${PREFIX}/lib/lib${LIB}.a
	rm ${PREFIX}/${INCLUDEDIR}/${LIB}.hpp
	rm ${PREFIX}/${INCLUDEDIR}/${LIB}_exception.hpp
	rm ${PREFIX}/lib/lib${LIB}.so.0.0.0
	rm ${PREFIX}/lib/lib${LIB}.so.0

