# Copyright 2012 Christian Eder
# 
# This file is part of F4RT, licensed under the GNU General Public
# License version 3. See COPYING for more information.

INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/include -I$(top_builddir)

#AM_CXXFLAGS = -g -Wall -DDEBUG -O0 -fomit-frame-pointer

SUBDIRS=

BASIC_HDRS =									\
	../include/f4rt-config.h \
	intrinsics.h

SOURCES=							\
	intrinsics.c

pkginclude_HEADERS =	\
	$(BASIC_HDRS)

libf4rt_la_CPPFLAGS= $(AM_CXXFLAGS) -c -Wall -pthread -D__F4RT_NDEBUG -DNDEBUG $(INCLUDES)

libf4rt_la_SOURCES=	\
	$(SOURCES)

libf4rt_la_LIBADD=
if ENABLE_DEBUG
AM_CPPFLAGS= -g3 -ggdb -Wall -pthread -D__F4RT_DEBUG -DDEBUG $(INCLUDES)
lib_LTLIBRARIES=libf4rt.la libf4rt_debug.la
libf4rt_debug_la_CPPFLAGS=$(AM_CPPFLAGS)
libf4rt_debug_la_SOURCES=	\
	$(SOURCES)
libf4rt_debug_la_LIBADD=
else
lib_LTLIBRARIES=libf4rt.la
endif

bin_PROGRAMS = f4rt

f4rt_CPPFLAGS= $(AM_CPPFLAGS)
f4rt_SOURCES= intrinsics.c intrinsics.h
#f4rt_LDADD=libf4rt.la
