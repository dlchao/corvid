### make changes accordingly ###
CC       = gcc
CPP      = g++
CLINKER  = gcc
CCLINKER = g++
MAKE     = make --no-print-directory
SHELL    = /bin/sh
CFLAGS   = -Wall -pedantic -DDSFMT_MEXP=19937 -msse2 -DHAVE_SSE2 
OPTI     = -O1
#OPTI = -pg # for profiling
LDFLAGS	= -lm
INCLUDES	= 
LIBS	= 
OBJS	= corvid.o epimodel.o params.o epimodelparameters.o dSFMT19937.o bnldev.o
DEFINES = -DVERBOSE

default: corvid

corvid: $(OBJS) Makefile
	$(CCLINKER) -o corvid $(OBJS) $(OPTI) $(LDFLAGS) $(LIBS) $(DEFINES)

R0corvid: $(OBJS) Makefile R0model.o R0model.h
	$(CCLINKER) -o R0corvid R0model.o epimodel.o params.o epimodelparameters.o dSFMT19937.o bnldev.o $(LDFLAGS) $(LIBS) $(DEFINES)

%.o: %.cpp epimodel.h epimodelparameters.h params.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c $<

dSFMT19937.o: dSFMT.c dSFMT.h dSFMT-params19937.h Makefile
	$(CC) $(CFLAGS) $(OPTI)  -std=c99 --param max-inline-insns-single=1800 -fno-strict-aliasing -Wmissing-prototypes -msse2 -DHAVE_SSE2 -DNDEBUG $(INCLUDES) $(DEFINES) -c dSFMT.c -o dSFMT19937.o

bnldev.o: bnldev.c bnldev.h Makefile
	$(CC) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c bnldev.c -o bnldev.o

zip: *.c *.cpp *.h Makefile
	cd ../..; zip corvid/corvid.zip corvid/README corvid/LICENSE corvid/gpl.txt corvid/HISTORY corvid/code/Makefile corvid/code/*.cpp corvid/code/*.c corvid/code/*.h corvid/corviddata/one-*dat corvid/corviddata/seattle-*dat corvid/corviddata/la-*dat corvid/configfiles/*

# makes directory ../corviddemo/ and runs several Seattle scenarios
demo: corvid
	rm -rf ../corviddemo
	mkdir ../corviddemo
	cd ../corviddemo; ln -s ../code/corvid .; ln -s ../corviddata/seattle* .; ln -s ../configfiles/config-seattle26* .
	cd ../corviddemo; ./corvid config-seattle26; ./corvid config-seattle26-closeallschools; ./corvid config-seattle26-isolation;  ./corvid config-seattle26-liberalleave; ./corvid config-seattle26-quarantine; ./corvid config-seattle26-workfromhome

# using emacs makes me feel old
emacs:
	emacs Makefile *.h *.c *.cpp &

clean:
	rm -f *.o corvid R0corvid *~
	rm -f Summary? Tracts? Log? Individuals?
	echo 0 > run-number
