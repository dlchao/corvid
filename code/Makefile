### make changes accordingly ###
CC       = gcc
CPP      = g++
CLINKER  = gcc
CCLINKER = g++
MAKE     = make --no-print-directory
SHELL    = /bin/sh
CFLAGS	 = -Wall -Werror -pedantic -DDSFMT_MEXP=19937 -msse2 -DHAVE_SSE2
# Check if OPTI not defined as an environment variable set the default
ifndef OPTI
OPTI     := -O3 -march=native
#OPTI = -pg # for profiling
endif
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
	cd ../..; zip corvid/corvid.zip corvid/README corvid/LICENSE corvid/gpl.txt corvid/HISTORY corvid/code/Makefile corvid/code/*.cpp corvid/code/*.c corvid/code/*.h corvid/corviddata/one-*dat corvid/corviddata/seattle-*dat corvid/corviddata/la-*dat corvid/sampleconfigfiles/*

# makes directory ../corvid-guide/ and runs several Seattle scenarios
guide: corvid
#	rm -rf ../corvid-guide
	mkdir -p ../corvid-guide
	cd ../corvid-guide; ln -s ../code/corvid .; ln -s ../corviddata/seattle* .; cp ../scripts/corvid-guide/* .
	cd ../corvid-guide; ./corvid config-seattle26; ./corvid config-seattle26-closeallschools; ./corvid config-seattle26-isolation;  ./corvid config-seattle26-liberalleave; ./corvid config-seattle26-quarantine; ./corvid config-seattle26-workfromhome; ./corvid config-seattle21-threshold;  ./corvid config-seattle26-threshold

# makes directory ../cdcmarch2020/ and runs several Seattle scenarios in response to the CDC coronavirus Modeling Team
# The runs will take hours. You might want to customize the run script to take advantage of computer resources better than my own.
cdcmarch2020: corvid
#	rm -rf ../cdcmarch2020 # you can erase this directory if you want to
	mkdir -p ../cdcmarch2020
	cd ../cdcmarch2020; ln -s ../code/corvid .; ln -s ../corviddata/seattle* .; cp ../scripts/cdcmodelingteam/run-cdc-2020-03-13.sh .; cp ../scripts/cdcmodelingteam/template*cdc-2020-03-13 .; cp ../scripts/cdcmodelingteam/config*cdc-2020-03-13 .; cp ../scripts/cdcmodelingteam/corvid-cdc-2020-03-13.Rmd .; chmod +x run-*.sh
	cd ../cdcmarch2020; ./run-cdc-2020-03-13.sh # run all the simulations for hours

# makes directory ../schoolclosuremarch2020/ and runs school closure
# The runs will take hours. You might want to customize the run script to take advantage of computer resources better than my own.
schoolclosuremarch2020: corvid Makefile
#	rm -rf ../schoolclosuremarch2020 # you can erase this directory if you want to
	mkdir -p ../schoolclosuremarch2020
	cd ../schoolclosuremarch2020; ln -s ../code/corvid .; ln -s ../corviddata/seattle* .; cp ../scripts/schoolclosure/* .; chmod +x run-*.sh
	cd ../schoolclosuremarch2020; ./run-schoolclosures-2020-03.sh # run all the simulations for hours

# makes directory ../layeredapril2020/ and runs layered NPI simulations
# The runs will take many hours. You might want to customize the run script to take advantage of computer resources better than my own.
layerednpisapril2020: corvid Makefile
#	rm -rf ../layerednpisapril2020 # you can erase this directory if you want to
	mkdir -p ../layerednpisapril2020
	cd ../layerednpisapril2020; ln -s ../code/corvid .; ln -s ../corviddata/seattle* .; cp ../scripts/layered/* .; chmod +x run-*.sh
	cd ../layerednpisapril2020; ./run-fig10.sh; ./run-fig1-top.sh; ./run-fig3-top.sh; ./run-fig5.sh; ./run-fig1-bottom.sh; ./run-fig3-bottom.sh; ./run-fig4.sh; mkdir -p results; mv Summary-* Log-* Tracts-* Individuals-* results # run all the simulations for hours

# using emacs makes me feel old
emacs:
	emacs Makefile *.h *.c *.cpp &

# does not clear out the guide or other directories
clean:
	rm -f *.o corvid R0corvid *~
	rm -f Summary? Tracts? Log? Individuals?

