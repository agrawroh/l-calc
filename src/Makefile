# Makefile for the L-function calculator version 1



#=============== c compiler ==================================
CC = g++

#cc = /Users/mrubinst/local/gcc4.4/bin/gcc
#CC = /Users/mrubinst/local/gcc4.4/bin/g++
#LD = /Users/mrubinst/local/gcc4.4/bin/g++



#=============== install directory (for 'make install') ======
INSTALL_DIR= /usr/local

#=============== location of various optional packages =======
# If you do not wish to use any of these packages, there is no need
# to change any of the lines below. By default they are not used.
# However, if you wish to use them (ex pari for some initializing an
# elliptic curve L-function) and do not have them installed
# in the standard location (/usr/local/include and /usr/local/lib)
# you can modify these variables.

#location of Bailey's include/qd/qd_real.h and include/qd/dd_real.h and
#of his lib/libqd.a
BAILEY_DIR = /usr/local
#BAILEY_dir=/Users/mrubinst/local/qd-2.3.8/install_dir
#on riemann.math
#/opt/qd_gnu

#location of PARI files: i.e. directory that has include/pari and lib with
PARI_DIR=/usr/local
#PARI_DIR=/Users/mrubinst/local/pari-2.3.4/install_dir
#on riemann.math
#/opt/pari_gnu

#location of libgmpfrxx libmpfr libgmp libgmpxx
#(for simplicity, these should all be in the same directory)
MPFR_LIBRARIES = /usr/local/lib

#===================== openmp ================================
#by default, openmp is turned off
#OPENMP_FLAG = -fopenmp

#=================== include files for L package =============
INCLUDEFILES:= -I../include

#================ operating system name ======================
OS_NAME := $(shell uname)

#===================== Choice of Doubles =====================
# Select what to use for doubles, and set appropriate variables.
# For example, typing 'make DOUBLE=dd', uses Bailey's quad doubles.
# Alternatively, setting DOUBLE=dd in the line below accomplishes
# the same.

# default (blank) is machine doubles. Other possibilities are
# Double = dd
# Double = qd
# Double = mpfr
# Double = ld
# These stand for Bailey's double double, Bailey's quad double, mpfr, and machine long doubles
# respectively.

DOUBLE =
LDFLAGS =

# 1) default is to use machine doubles
ifeq ($(DOUBLE),)
    PREPROCESSOR_DEFINE = -DUSE_DOUBLE
    EXTENSION =
    EXTRA := --fast-math -fno-finite-math-only
endif

# 2) Bailey's double double. Gives about 30 digits precision.
# Is about ten times faster than the c++ version of mpfr that I am using (gmpfrxx).
ifeq ($(DOUBLE),dd)
    PREPROCESSOR_DEFINE = -DUSE_BAILEY_DD
    INCLUDEFILES += -I$(BAILEY_DIR)/include
    LDFLAGS += -L$(BAILEY_DIR)/lib -lqd
    EXTENSION = _dd
    EXTRA:=
endif

# 3) Bailey's quad double. Gives around 60 Digits precision.
# Is still faster than c++ mpfr
ifeq ($(DOUBLE),qd)
    PREPROCESSOR_DEFINE = -DUSE_BAILEY_QD
    INCLUDEFILES += -I$(BAILEY_DIR)/include
    LDFLAGS += -L$(BAILEY_DIR)/lib -lqd
    EXTENSION = _qd
    EXTRA:=
endif

# 4) mpfr package. Should be used if more than 60 digits accuracy is needed.
ifeq ($(DOUBLE),mpfr)
    PREPROCESSOR_DEFINE = -DUSE_MPFR
    LDFLAGS += -L$(MPFR_LIBRARIES) -lgmpfrxx -lmpfr -lgmp -lgmpxx -lm
    #EXTRA_FLAGS= gmpfrxx.o mpfr_mul_d.o -lmpfr -lgmp -lgmpxx -lm
    EXTENSION = _mpfr
    EXTRA:=
endif

# 5) mpfrc++ package. Should be used if more than 60 digits accuracy is needed.
ifeq ($(DOUBLE),mpfrcpp)
    PREPROCESSOR_DEFINE = -DUSE_MPFRCPP
    LDFLAGS += mpreal.o dlmalloc.o -L$(MPFR_LIBRARIES) -lmpfr -lgmp -lm
    #g++ example.o mpreal.o dlmalloc.o /usr/local/lib/libmpfr.a /usr/local/lib/libgmp.a -o example
    EXTENSION = _mpfrc++
    EXTRA:=
endif

# 6) machine long doubles
ifeq ($(DOUBLE),ld)
    PREPROCESSOR_DEFINE = -DUSE_LONG_DOUBLE
    EXTENSION = _ld
    EXTRA := --fast-math
endif


#================ pari ======================================
# Comment out the following line to remove the use of pari's
# elliptic curve routines. Doing so disables lcalc's -e commandline option.
# For the -e lcalc option to work you need pari installed on your machine.
PARI_DEFINE = -DINCLUDE_PARI

ifeq ($(PARI_DEFINE),-DINCLUDE_PARI)
    LDFLAGS += -L$(PARI_DIR)/lib -lpari
endif

#==================== misc extra compiling options ============

#EXTRA += -ftree-vectorize  -ftree-vectorizer-verbose=5
#-g for Apple's shark program to show me source code in profiling
EXTRA += -g
#EXTRA += -fomit-frame-pointer -mtune=core2 #-march=native

MACHINE_SPECIFIC_FLAGS :=
#warning- fast math causes trouble with Bailey's double/quad double package and also gmp:
#MACHINE_SPECIFIC_FLAGS += -ffast-math
#on my old G5:
#MACHINE_SPECIFIC_FLAGS +=  -mpowerpc -mpowerpc64 -m64

#================== disable a few of -Wall warnings ==========
DISABLE_WARNINGS:= -Wno-unused-parameter -Wno-uninitialized

ifneq ($(OPENMP_FLAG),-fopenmp)
    DISABLE_WARNINGS+= -Wno-unknown-pragmas
endif

#==================== set the compiling flags ============
#CCFLAGS = -fPIC -Wall $(DISABLE_WARNINGS) $(OPENMP_FLAG) $(PREPROCESSOR_DEFINE) $(MACHINE_SPECIFIC_FLAGS) $(EXTRA)
CCFLAGS = -O3 -fPIC -Wall $(DISABLE_WARNINGS) $(OPENMP_FLAG) $(PREPROCESSOR_DEFINE) $(MACHINE_SPECIFIC_FLAGS) $(EXTRA)

#CCFLAGS =  -Wa,-W -O2 $(OPENMP_FLAG) -Wno-deprecated $(PREPROCESSOR_DEFINE) $(MACHINE_SPECIFIC_FLAGS) $(EXTRA)
#CCFLAGS =  -Wall -O2 $(OPENMP_FLAG) $(PREPROCESSOR_DEFINE) $(MACHINE_SPECIFIC_FLAGS) $(EXTRA)


#================ various linking flags =======================
#For mac os x shared libraries are hanled with the -dynamiclib option
ifeq ($(OS_NAME),Darwin)
    DYN_OPTION=dynamiclib
else
    DYN_OPTION=shared

    LDFLAGS += -Xlinker -rpath -Xlinker . #for linking to Liblfunction within the current directory

    ifeq ($(PARI_DEFINE),-DINCLUDE_PARI)
       LDFLAGS += -Xlinker -export-dynamic #not sure why pari calls these when linking but on the web I found
       #'Libtool provides the `-export-dynamic' link flag (see section Link mode), which does this declaration.
       #You need to use this flag if you are linking a shared library that will be dlopened'
       LDFLAGS += -Xlinker -rpath -Xlinker $(PARI_DIR)/lib
    endif
endif


ifeq ($(OPENMP_FLAG),-fopenmp)
    LDFLAGS += -lgomp
endif

#Xlinker NOTES:
#At runtime, the computer attempts to load the pari shared library, and if it isn't in a standard
#location, we can do two things.
#1) One page suggested:
#    -Xlinker -rpath -Xlinker /usr/local/pari/pari-2.1.5/lib
#2) The other option, not recommended, is to type at the unix prompt:
#    LD_LIBRARY_PATH=/usr/local/pari/pari-2.1.5/lib:$LD_LIBRARY_PATH
#    export LD_LIBRARY_PATH
#If this is not done correctly, at runntime one gets the error
#    ./lcalc: error while loading shared libraries: libpari.so.1: cannot open shared
#    object file: No such file or directory
#One can list, after compiling, dynamic dependencies with the command: ldd lcalc and it will
#become clear which libraries the computer can find.

#LDFLAGS += -L/Users/mrubinst/local/gcc4.4/lib

#===================== object files ==============================
#object files for the libLfunction library
#OBJ_L = Lglobals.o Lmisc.o Lgamma.o Lriemannsiegel.o Lriemannsiegel_blfi.o Lnumbertheory.o Lelliptic.o Ldokchitser.o
OBJ_L = Lglobals.o Lmisc.o Lgamma.o Lriemannsiegel.o Lnumbertheory.o Lelliptic.o Ldokchitser.o

#object files for the command line program
#OBJECTS := $(OBJ_L) Lcommandline_globals.o Lcommandline_misc.o Lcommandline_values_zeros.o
OBJECTS := Lcommandline_globals.o Lcommandline_misc.o Lcommandline_values_zeros.o
OBJECTS += Lcommandline_twist.o Lcommandline.o cmdline.o

#================ source files (used by 'make depend') ===========
#SRCS = Lcommandline.cc Lcommandline_globals.cc Lcommandline_misc.cc Lcommandline_twist.cc Lcommandline_values_zeros.cc Lgamma.cc Lglobals.cc Lmisc.cc Lriemannsiegel.cc Lriemannsiegel_blfi.cc Lnumbertheory.cc Lelliptic.cc cmdline.c
SRCS = Lcommandline.cc Lcommandline_globals.cc Lcommandline_misc.cc Lcommandline_twist.cc Lcommandline_values_zeros.cc Lgamma.cc Lglobals.cc Lmisc.cc Lriemannsiegel.cc Lnumbertheory.cc Lelliptic.cc cmdline.c


#========================= make ==================================
all:
	$(MAKE) library
	$(MAKE) lcalc
	$(MAKE) examples
#	$(MAKE) sacha
	$(MAKE) test
#	make library
#	make lcalc
#	make examples
#	make sacha
#	make test
#	make find_L
#	make test_moment
#	make print_vars


library: $(OBJ_L)
	$(CC) -$(DYN_OPTION) $(LDFLAGS) -o libLfunction$(EXTENSION).so $(OBJ_L)
#	$(CC) -$(DYN_OPTION) -o libLfunction$(EXTENSION).so $(OBJ_L)

lcalc: $(OBJECTS)
	$(CC) $(CCFLAGS) $(OBJECTS) -L. -lLfunction$(EXTENSION) $(LDFLAGS) -o lcalc$(EXTENSION)
#	$(CC) $(CCFLAGS) $(INCLUDEFILES) $(OBJECTS) $(LDFLAGS) -o lcalc$(EXTENSION)

Lglobals.o: Lglobals.cc
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -I$(PARI_DIR)/include $(PARI_DEFINE) -c Lglobals.cc

Lelliptic.o: Lelliptic.cc
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -I$(PARI_DIR)/include $(PARI_DEFINE) -c Lelliptic.cc

Lcommandline_twist.o: Lcommandline_twist.cc
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -I$(PARI_DIR)/include $(PARI_DEFINE) -c Lcommandline_twist.cc

Lcommandline.o: Lcommandline.cc
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -I$(PARI_DIR)/include $(PARI_DEFINE) -c Lcommandline.cc

.cc.o:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -c $<

.c.o:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -c $<
#$(CC) $(CCFLAGS) $(INCLUDEFILES) -DHAVE_LONG_LONG -c cmdline.c

examples:
#	$(CC) $(CCFLAGS) $(INCLUDEFILES) example_programs/example2.cc -L. -lLfunction$(EXTENSION) $(LDFLAGS) -o example_programs/example2
	$(CC) $(CCFLAGS) $(INCLUDEFILES) example_programs/example.cc -L. -lLfunction$(EXTENSION) $(LDFLAGS) -o example_programs/example

sacha:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -I$(PARI_DIR)/include $(PARI_DEFINE) sacha.cc -L. -lLfunction$(EXTENSION) $(LDFLAGS) -o sacha$(EXTENSION)


test:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) -I$(PARI_DIR)/include $(PARI_DEFINE) example_programs/test_accuracy.cc -L. -lLfunction$(EXTENSION) $(LDFLAGS) -o example_programs/test_accuracy$(EXTENSION)

proc:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) example_programs/proc.cc $(LDFLAGS) libLfunction$(EXTENSION).so -o example_programs/proc

test_moment:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) example_programs/test_moment.cc $(LDFLAGS) libLfunction$(EXTENSION).so -o example_programs/test

find_L:
	$(CC) $(CCFLAGS) $(INCLUDEFILES) find_L_functions/find_L_functions.cc $(LDFLAGS) libLfunction$(EXTENSION).so -o find_L_functions/find_L

clean:
	rm -f L*.o cmdline.o example_programs/example example_programs/test_accuracy sacha
#	rm -f L*.o cmdline.o lcalc$(EXTENSION) example_programs/example example_programs/test_accuracy example_programs/sacha

install:
	cp -f lcalc$(EXTENSION) $(INSTALL_DIR)/bin/.
	cp -f libLfunction$(EXTENSION).so $(INSTALL_DIR)/lib/.
	cp -rf ../include $(INSTALL_DIR)/include/Lfunction

depend:
	makedepend -- $(CCFLAGS) -Y../include -- $(SRCS)

print_vars:
	@echo OS_NAME = $(OS_NAME)
	@echo LDFLAGS = $(LDFLAGS)


# DO NOT DELETE

Lcommandline.o: ../include/Lcommandline.h ../include/L.h
Lcommandline.o: ../include/Lglobals.h ../include/Lcomplex.h
Lcommandline.o: ../include/Lcommon.h ../include/Lmisc.h ../include/Lgamma.h
Lcommandline.o: ../include/Lriemannsiegel.h ../include/Lriemannsiegel_blfi.h
Lcommandline.o: ../include/Lnumbertheory.h ../include/Lelliptic.h
Lcommandline.o: ../include/Ldirichlet_series.h ../include/Lprint.h
Lcommandline.o: ../include/Lnumberzeros.h ../include/Lgram.h
Lcommandline.o: ../include/Lvalue.h ../include/Lfind_zeros.h
Lcommandline.o: ../include/Ldokchitser.h ../include/Lexplicit_formula.h
Lcommandline.o: ../include/Lcommandline_globals.h
Lcommandline.o: ../include/Lcommandline_misc.h
Lcommandline.o: ../include/Lcommandline_twist.h
Lcommandline.o: ../include/Lcommandline_values_zeros.h ../include/cmdline.h
Lcommandline_globals.o: ../include/Lcommandline_globals.h ../include/L.h
Lcommandline_globals.o: ../include/Lglobals.h ../include/Lcomplex.h
Lcommandline_globals.o: ../include/Lcommon.h ../include/Lmisc.h
Lcommandline_globals.o: ../include/Lgamma.h ../include/Lriemannsiegel.h
Lcommandline_globals.o: ../include/Lriemannsiegel_blfi.h
Lcommandline_globals.o: ../include/Lnumbertheory.h ../include/Lelliptic.h
Lcommandline_globals.o: ../include/Ldirichlet_series.h ../include/Lprint.h
Lcommandline_globals.o: ../include/Lnumberzeros.h ../include/Lgram.h
Lcommandline_globals.o: ../include/Lvalue.h ../include/Lfind_zeros.h
Lcommandline_globals.o: ../include/Ldokchitser.h
Lcommandline_globals.o: ../include/Lexplicit_formula.h
Lcommandline_misc.o: ../include/Lcommandline_misc.h ../include/L.h
Lcommandline_misc.o: ../include/Lglobals.h ../include/Lcomplex.h
Lcommandline_misc.o: ../include/Lcommon.h ../include/Lmisc.h
Lcommandline_misc.o: ../include/Lgamma.h ../include/Lriemannsiegel.h
Lcommandline_misc.o: ../include/Lriemannsiegel_blfi.h
Lcommandline_misc.o: ../include/Lnumbertheory.h ../include/Lelliptic.h
Lcommandline_misc.o: ../include/Ldirichlet_series.h ../include/Lprint.h
Lcommandline_misc.o: ../include/Lnumberzeros.h ../include/Lgram.h
Lcommandline_misc.o: ../include/Lvalue.h ../include/Lfind_zeros.h
Lcommandline_misc.o: ../include/Ldokchitser.h ../include/Lexplicit_formula.h
Lcommandline_misc.o: ../include/Lcommandline_globals.h
Lcommandline_twist.o: ../include/Lcommandline_twist.h ../include/L.h
Lcommandline_twist.o: ../include/Lglobals.h ../include/Lcomplex.h
Lcommandline_twist.o: ../include/Lcommon.h ../include/Lmisc.h
Lcommandline_twist.o: ../include/Lgamma.h ../include/Lriemannsiegel.h
Lcommandline_twist.o: ../include/Lriemannsiegel_blfi.h
Lcommandline_twist.o: ../include/Lnumbertheory.h ../include/Lelliptic.h
Lcommandline_twist.o: ../include/Ldirichlet_series.h ../include/Lprint.h
Lcommandline_twist.o: ../include/Lnumberzeros.h ../include/Lgram.h
Lcommandline_twist.o: ../include/Lvalue.h ../include/Lfind_zeros.h
Lcommandline_twist.o: ../include/Ldokchitser.h ../include/Lexplicit_formula.h
Lcommandline_twist.o: ../include/Lcommandline_globals.h
Lcommandline_values_zeros.o: ../include/L.h ../include/Lglobals.h
Lcommandline_values_zeros.o: ../include/Lcomplex.h ../include/Lcommon.h
Lcommandline_values_zeros.o: ../include/Lmisc.h ../include/Lgamma.h
Lcommandline_values_zeros.o: ../include/Lriemannsiegel.h
Lcommandline_values_zeros.o: ../include/Lriemannsiegel_blfi.h
Lcommandline_values_zeros.o: ../include/Lnumbertheory.h
Lcommandline_values_zeros.o: ../include/Lelliptic.h
Lcommandline_values_zeros.o: ../include/Ldirichlet_series.h
Lcommandline_values_zeros.o: ../include/Lprint.h ../include/Lnumberzeros.h
Lcommandline_values_zeros.o: ../include/Lgram.h ../include/Lvalue.h
Lcommandline_values_zeros.o: ../include/Lfind_zeros.h
Lcommandline_values_zeros.o: ../include/Ldokchitser.h
Lcommandline_values_zeros.o: ../include/Lexplicit_formula.h
Lcommandline_values_zeros.o: ../include/Lcommandline_values_zeros.h
Lcommandline_values_zeros.o: ../include/Lcommandline_globals.h
Lgamma.o: ../include/Lgamma.h ../include/Lglobals.h ../include/Lcomplex.h
Lgamma.o: ../include/Lcommon.h ../include/Lmisc.h
Lglobals.o: ../include/Lglobals.h ../include/Lcomplex.h ../include/Lcommon.h
Lmisc.o: ../include/Lmisc.h ../include/Lglobals.h ../include/Lcomplex.h
Lmisc.o: ../include/Lcommon.h
Lriemannsiegel.o: ../include/L.h ../include/Lglobals.h ../include/Lcomplex.h
Lriemannsiegel.o: ../include/Lcommon.h ../include/Lmisc.h ../include/Lgamma.h
Lriemannsiegel.o: ../include/Lriemannsiegel.h
Lriemannsiegel.o: ../include/Lriemannsiegel_blfi.h ../include/Lnumbertheory.h
Lriemannsiegel.o: ../include/Lelliptic.h ../include/Ldirichlet_series.h
Lriemannsiegel.o: ../include/Lprint.h ../include/Lnumberzeros.h
Lriemannsiegel.o: ../include/Lgram.h ../include/Lvalue.h
Lriemannsiegel.o: ../include/Lfind_zeros.h ../include/Ldokchitser.h
Lriemannsiegel.o: ../include/Lexplicit_formula.h
#Lriemannsiegel_blfi.o: ../include/L.h ../include/Lglobals.h
#Lriemannsiegel_blfi.o: ../include/Lcomplex.h ../include/Lcommon.h
#Lriemannsiegel_blfi.o: ../include/Lmisc.h ../include/Lgamma.h
#Lriemannsiegel_blfi.o: ../include/Lriemannsiegel.h
#Lriemannsiegel_blfi.o: ../include/Lriemannsiegel_blfi.h
#Lriemannsiegel_blfi.o: ../include/Lnumbertheory.h ../include/Lelliptic.h
#Lriemannsiegel_blfi.o: ../include/Ldirichlet_series.h ../include/Lprint.h
#Lriemannsiegel_blfi.o: ../include/Lnumberzeros.h ../include/Lgram.h
#Lriemannsiegel_blfi.o: ../include/Lvalue.h ../include/Lfind_zeros.h
#Lriemannsiegel_blfi.o: ../include/Ldokchitser.h
#Lriemannsiegel_blfi.o: ../include/Lexplicit_formula.h
Lnumbertheory.o: ../include/Lnumbertheory.h ../include/Lglobals.h
Lnumbertheory.o: ../include/Lcomplex.h ../include/Lcommon.h
Lelliptic.o: ../include/Lelliptic.h ../include/L.h ../include/Lglobals.h
Lelliptic.o: ../include/Lcomplex.h ../include/Lcommon.h ../include/Lmisc.h
Lelliptic.o: ../include/Lgamma.h ../include/Lriemannsiegel.h
Lelliptic.o: ../include/Lriemannsiegel_blfi.h ../include/Lnumbertheory.h
Lelliptic.o: ../include/Ldirichlet_series.h ../include/Lprint.h
Lelliptic.o: ../include/Lnumberzeros.h ../include/Lgram.h ../include/Lvalue.h
Lelliptic.o: ../include/Lfind_zeros.h ../include/Ldokchitser.h
Lelliptic.o: ../include/Lexplicit_formula.h
cmdline.o: ../include/getopt.h ../include/cmdline.h
