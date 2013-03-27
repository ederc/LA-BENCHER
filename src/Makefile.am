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
	f4rt-kaapi.h							\
	f4rt-omp.h								\
	f4rt-pthread.h						\
	f4rt-seq.h								\
	f4rt-tbb.h								\
	matrix.h									\
	dense-mat-mult.h

SOURCES=										\
	intrinsics.c							\
	f4rt.cpp									\
	f4rt-kaapi.cpp						\
	f4rt-omp.cpp							\
	f4rt-pthread.cpp					\
	f4rt-seq.cpp							\
	f4rt-tbb.cpp							\
	matrix.cpp								\
	dense-mat-mult.cpp

pkginclude_HEADERS =	\
	$(BASIC_HDRS)

libf4rt_la_CPPFLAGS= $(AM_CXXFLAGS) -D__F4RT_NDEBUG -DNDEBUG $(INCLUDES)
libf4rt_la_CXXFLAGS= $(AM_CXXFLAGS) -std=gnu++0x -Wall -D__F4RT_NDEBUG -DNDEBUG $(INCLUDES)

libf4rt_la_SOURCES=	\
	$(SOURCES)

libf4rt_la_LIBADD=
if ENABLE_DEBUG
lib_LTLIBRARIES=libf4rt.la libf4rt_debug.la
libf4rt_debug_la_CPPFLAGS=$(AM_CXXFLAGS) -g3 -ggdb -Wall -D__F4RT_DEBUG -DDEBUG $(INCLUDES)
libf4rt_debug_la_CXXFLAGS=$(AM_CXXFLAGS) -g3 -ggdb -Wall -std=gnu++0x -D__F4RT_DEBUG -DDEBUG $(INCLUDES)
libf4rt_debug_la_SOURCES=	\
	$(SOURCES)
libf4rt_debug_la_LIBADD=
else
lib_LTLIBRARIES=libf4rt.la
endif

bin_PROGRAMS = intrinsics f4rt

intrinsics_CPPFLAGS= $(AM_CXXFLAGS)
intrinsics_SOURCES= intrinsics.c intrinsics.h
intrinsics_LDADD=libf4rt.la

f4rt_CXXFLAGS= $(AM_CXXFLAGS) -std=gnu++0x
f4rt_SOURCES=  				\
	f4rt.cpp 						\
	matrix.cpp 					\
	f4rt-kaapi.cpp 			\
	f4rt-omp.cpp 				\
	f4rt-pthread.cpp 		\
	f4rt-seq.cpp				\
	f4rt-tbb.cpp				\
	dense-mat-mult.cpp
