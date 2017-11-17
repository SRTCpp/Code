PROGRAM = $@
JDIR=~/working/c++/Jcube
NR = ~/working/nr/
BINDIR = ~/bin
LIBS =  -lcspice -lcsupport -lm -lboost_program_options -lcommandl \
 -lcfitsio -lnr -lgdal\
 -lc++ -lJcube `Magick++-config --ldflags --libs` 
RM = /bin/rm -f
CCver = 49
CC = g++${CCver} -std=c++11 -fopenmp -g -I /usr/local/include -Wno-deprecated `Magick++-config --cxxflags --cppflags` 
LIBLOC=~/lib
LFLAGS= -Wl,-rpath,/usr/local/lib/gcc${CCver} -pthread
SPICELOC=/usr/local
TOPOSPICELOC=~/working/c++/Jcube/titantopo/spice
TOPOLOC=~/working/c++/Jcube/titantopo
RADTRANO= SRTC++.o photons.o atmolayers.o atmospheres.o geomvector.o \
   photongenerators.o \
   phasefunctions.o atmozones.o photontraverse.o detectors.o 
RADTRANH= SRTC++.h photons.h atmolayers.h atmospheres.h geomvector.h \
   photongenerators.h \
   phasefunctions.h atmozones.h photontraverse.h detectors.h \
	test_GrahamCompare.h SRTC++_testsuite.h
TESTSUITEO= test_LambertSurface.o test_ChandraAtmo.o \
  chandra_test.o test_SebastienCompare.o \
  test_GrahamCompare.o test_DiffusiveReflectance.o diffusive_reflectance.o 
BENCHO= 
TESTSUITEH= STRC++_testsuite.h chandra_test.h  
LIBDIRS= ${LFLAGS} -L ~/lib -L/usr/local/lib/gcc${CCver} -L/usr/local/lib  \
 -L/usr/local/lib/mysql -L ${LIBLOC} -L/usr/lib

.SUFFIXES  : .c++ .o
.c++.o:
	${CC} -c $<
%.d: %.c++
	$(CC) -MM -MD $<
							 
include $(RADTRANO:.o=.d)
include $(TESTSUITEO:.o=.d)

clean:
	${RM} *.o *core ${PROGRAM}*.o

#${PROGRAM}.o: $(PROGRAM).c++ ${JDIR}/Jcube.h ${JDIR}/Junit.h \
 ${JDIR}/Jvalue.h ${JDIR}/Jdiffeq.h ${JDIR}/Jdiffeq.h Jcube_lib \
 ${LIBLOC}/libJcube.a ${LIBLOC}/libnr.a $(RADTRANH) $(RADTRANO) $@.d
#	$(CC) -c $(PROGRAM).c++ -o ${PROGRAM}.o

#$@: ${PROGRAM}.o Jcube_lib ${LIBLOC}/libJcube.a ${LIBLOC}/libnr.a $(RADTRANO)
#	$(CC) $< ../JIM.o $(RADTRANO) $(TESTSUITEO) \
	$(LIBDIRS) -o $(BINDIR)/${PROGRAM} \
 	$(LIBS) ${LIBS}

diskintegrate: diskintegrate.o ${LIBLOC}/libJcube.a ${LIBLOC}/libnr.a $(RADTRANO)
	$(CC) $< ../JIM.o $(RADTRANO) \
	$(LIBDIRS) -o $(BINDIR)/$@ \
 	$(LIBS) ${LIBS}
	
SRTC++test: SRTC++test.o ${LIBLOC}/libJcube.a ${LIBLOC}/libnr.a $(RADTRANO)
	$(CC) $< ../JIM.o $(RADTRANO) \
	$(LIBDIRS) -o $(BINDIR)/$@ \
 	$(LIBS) ${LIBS}
		
SRTC++_testsuite: SRTC++_testsuite.o ${LIBLOC}/libJcube.a ${LIBLOC}/libnr.a $(RADTRANO) $(TESTSUITEO) 
	$(CC) $< ../JIM.o $(RADTRANO) $(TESTSUITEO) \
	$(LIBDIRS) -o $(BINDIR)/$@ \
 	$(LIBS) ${LIBS}
		
SRTC++_benchmark: SRTC++_benchmark.o ${LIBLOC}/libJcube.a ${LIBLOC}/libnr.a $(RADTRANO) $(BENCHO) 
	$(CC) $< ../JIM.o $(RADTRANO) $(BENCHO) \
	$(LIBDIRS) -o $(BINDIR)/$@ \
 	$(LIBS) ${LIBS}
		
syntheticimageintegrate: syntheticimageintegrate.o ${LIBLOC}/libJcube.a ${LIBLOC}/libnr.a $(RADTRANO)
	$(CC) $< ../JIM.o $(RADTRANO) \
	$(LIBDIRS) -o $(BINDIR)/$@ \
 	$(LIBS) ${LIBS}

