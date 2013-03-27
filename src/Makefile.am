# Copyright 2012 Christian Eder
# 
# This file is part of F4RT, licensed under the GNU General Public
# License version 3. See COPYING for more information.

AUTOMAKE_OPTIONS = subdir-objects

INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/include -I$(top_builddir)

#AM_CXXFLAGS = -g -Wall -DDEBUG -O0 -fomit-frame-pointer


SUBDIRS=

BASIC_HDRS =								\
	../include/f4rt-config.h 	\
	intrinsics.h							\
	kaapi/mat-mult-kaapi.h		\
	omp/mat-mult-omp.h				\
	pthrd/mat-mult-pthrd.h		\
	seq/mat-mult-seq.h				\
	tbb/mat-mult-tbb.h				\
	matrix.h									\
	mat-gen.h									\
	mat-mult.h

SOURCES=										\
	intrinsics.c							\
	f4rt.cpp									\
	kaapi/mat-mult-kaapi.cpp	\
	omp/mat-mult-omp.cpp			\
	pthrd/mat-mult-pthrd.cpp	\
	seq/mat-mult-seq.cpp			\
	tbb/mat-mult-tbb.cpp			\
	matrix.cpp								\
	mat-gen.cpp								\
	mat-mult.cpp

pkginclude_HEADERS =	\
	$(BASIC_HDRS)

libf4rt_la_CPPFLAGS	= $(AM_CXXFLAGS) -D__F4RT_NDEBUG -DNDEBUG $(INCLUDES)
libf4rt_la_CXXFLAGS	= $(AM_CXXFLAGS) -std=gnu++0x -Wall -D__F4RT_NDEBUG -DNDEBUG $(INCLUDES)

libf4rt_la_SOURCES	=	\
	$(SOURCES)

libf4rt_la_LIBADD=
if ENABLE_DEBUG
lib_LTLIBRARIES	=	libf4rt.la libf4rt_debug.la
libf4rt_debug_la_CPPFLAGS	=	$(AM_CXXFLAGS) -g3 -ggdb -Wall -D__F4RT_DEBUG -DDEBUG $(INCLUDES)
libf4rt_debug_la_CXXFLAGS	=	$(AM_CXXFLAGS) -g3 -ggdb -Wall -std=gnu++0x -D__F4RT_DEBUG -DDEBUG $(INCLUDES)
libf4rt_debug_la_SOURCES	=	\
	$(SOURCES)
libf4rt_debug_la_LIBADD=
else
lib_LTLIBRARIES	=	libf4rt.la
endif

bin_PROGRAMS = intrinsics f4rt

intrinsics_CPPFLAGS= $(AM_CXXFLAGS)
intrinsics_SOURCES= intrinsics.c intrinsics.h
intrinsics_LDADD=libf4rt.la

f4rt_CXXFLAGS= $(AM_CXXFLAGS) -std=gnu++0x
f4rt_SOURCES=  							\
	f4rt.cpp 									\
	matrix.cpp 								\
	kaapi/mat-mult-kaapi.cpp	\
	omp/mat-mult-omp.cpp			\
	pthrd/mat-mult-pthrd.cpp	\
	seq/mat-mult-seq.cpp			\
	tbb/mat-mult-tbb.cpp			\
	mat-gen.cpp								\
	mat-mult.cpp

CLEANFILES=				\
	../obj/*				\
	kaapi/.dirstamp	\
	pthrd/.dirstamp	\
	omp/.dirstamp		\
	tbb/.dirstamp		\
	seq/.dirstamp		

DISTCLEANFILES=		\
	../obj/*				\
	kaapi/.dirstamp	\
	pthrd/.dirstamp	\
	omp/.dirstamp		\
	tbb/.dirstamp		\
	seq/.dirstamp		

%/.dirstamp:
	mkdir -p $(@D)
	touch $@

all:
	mkdir -p ../obj
	mv *.o ../obj
	mv *.lo ../obj
	mv kaapi/*.o ../obj
	mv kaapi/*.lo ../obj
	mv omp/*.o ../obj
	mv omp/*.lo ../obj
	mv pthrd/*.o ../obj
	mv pthrd/*.lo ../obj
	mv seq/*.o ../obj
	mv seq/*.lo ../obj
	mv tbb/*.o ../obj
	mv tbb/*.lo ../obj
