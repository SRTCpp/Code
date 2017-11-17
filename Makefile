JDIR=~/working/c++/Jcube
NR = ~/working/nr/
JCUBEO = Jcube.o Jio.o Joperators.o Jmanip.o Jfunc.o Jmisc.o JIM.o Jnr.o Junit.o \
Jprocessors.o Jfits.o Jdate.o Jdraw.o Jalgorithms.o Jaddons.o Jassociated.o Jprojections.o Jangle.o
LIBS = -lm -lnr -lsqlplus -lCCfits -lcfitsio -lQt3Support -lQtSql -lQtXml -lQtNetwork -lQtGui -lQtCore -lgdal
RM = /bin/rm -f
CCver = 49
CC = g++${CCver} -fopenmp -std=c++11 -Wno-deprecated -g  ${INCLS} `Magick++-config --cxxflags --cppflags`
INCLS= -I /usr/local/include/mysql -I /usr/local/include -I/usr/local/share/qt/mkspecs/freebsd-g++ -I/usr/local/include/qt4/QtCore -I/usr/local/include/qt4/QtGui -I/usr/local/include/qt4/Qt3Support -I/usr/local/include/qt4
LIBLOC=~/lib
LFLAGS=-Wl,-rpath,/usr/local/lib/gcc${CCver}
LIBDIRS= -L /usr/local/lib/gcc${CCver} -L/usr/local/lib -L/usr/X11R6/lib -L/usr/local/lib/mysql -L ${LIBLOC} -L/usr/local/lib/qt4

.SUFFIXES  : .c++ .o
.c++.o:
	${CC} -c $<



all:	nr_lib ${LIBLOC}/libJcube.a

clean:
	${RM} *core *.o

nr_lib:
	cd ${NR} ; gmake

Jprocessors.o:	Jprocessors.c++ Jcube.h

Junit.o:	  		Junit.c++ Junit.h

Jnr.o:         Jnr.c++ Jnr.h Jcube.h

JIM.o:         JIM.c++ Jcube.h

Jfunc.o:       Jfunc.c++ Jcube.h Jangle.h

Jio.o:         Jio.c++ Jcube.h

Joperators.o:  Joperators.c++ Jcube.h

Jmanip.o:		Jmanip.c++ Jcube.h

Jmisc.o:			Jmisc.c++ Jcube.h

Jcube.o:			Jcube.c++ Jcube.h

Jfits.o:			Jfits.c++ Jcube.h

Jdate.o:			Jdate.c++ Jdate.h

Jdraw.o:			Jdraw.c++ Jcube.h

Jaddons.o:			Jaddons.c++ Jcube.h Jaddons.h

Jalgorithms.o:			Jalgorithms.c++ Jcube.h

Jassociated.o:			Jassociated.c++ Jcube.h Jassociated.h

Jprojections.o:		Jprojections.c++ Jcube.h Jprojections.h 

Jangle.o:		Jangle.c++ Jangle.h 

${LIBLOC}/libJcube.a:		${JCUBEO}
	ar r ${LIBLOC}/libJcube.a ${JCUBEO}
