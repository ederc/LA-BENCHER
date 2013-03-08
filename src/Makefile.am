# Copyright 2012 Christian Eder
# 
# This file is part of F4RT, licensed under the GNU General Public
# License version 3. See COPYING for more information.

INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/include -I$(top_builddir)

#AM_CXXFLAGS = -g -Wall -DDEBUG -O0 -fomit-frame-pointer

SUBDIRS=

BASIC_HDRS =								\
	../include/f4rt-config.h 	\
	intrinsics.h							\
	matrix-omp.h							\
	matrix-seq.h							\
	matrix-tbb.h							\
	matrix.h									

SOURCES=										\
	intrinsics.c							\
	dense-matrix-mult.cpp			\
	matrix-omp.cpp						\
	matrix-seq.cpp						\
	matrix-tbb.cpp						\
	matrix.cpp

pkginclude_HEADERS =	\
	$(BASIC_HDRS)

libf4rt_la_CPPFLAGS= $(AM_CPPFLAGS) -Wall -pthread -D__F4RT_NDEBUG -DNDEBUG $(INCLUDES)
libf4rt_la_CXXFLAGS= $(AM_CPPFLAGS) -std=gnu++0x -Wall -pthread -fopenmp -D__F4RT_NDEBUG -DNDEBUG $(INCLUDES)

libf4rt_la_SOURCES=	\
	$(SOURCES)

libf4rt_la_LIBADD=
if ENABLE_DEBUG
lib_LTLIBRARIES=libf4rt.la libf4rt_debug.la
libf4rt_debug_la_CPPFLAGS=$(AM_CPPFLAGS) -g3 -ggdb -Wall -pthread -D__F4RT_DEBUG -DDEBUG $(INCLUDES)
libf4rt_debug_la_CXXFLAGS=$(AM_CPPFLAGS) -g3 -ggdb -Wall -std=gnu++0x -pthread -fopenmp -D__F4RT_DEBUG -DDEBUG $(INCLUDES)
libf4rt_debug_la_SOURCES=	\
	$(SOURCES)
libf4rt_debug_la_LIBADD=
else
lib_LTLIBRARIES=libf4rt.la
endif

bin_PROGRAMS = f4rt dense-mult

f4rt_CPPFLAGS= $(AM_CPPFLAGS)
f4rt_SOURCES= intrinsics.c intrinsics.h
f4rt_LDADD=libf4rt.la

dense_mult_CXXFLAGS= $(AM_CPPFLAGS) -ltbb -std=gnu++0x -fopenmp -Wall 
dense_mult_SOURCES= dense-matrix-mult.cpp matrix.cpp matrix-tbb.cpp matrix-omp.cpp matrix-seq.cpp
