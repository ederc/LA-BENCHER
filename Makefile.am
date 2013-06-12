# Copyright 2012 Christian Eder
# 
# This file is part of LA-BENCHER, licensed under the GNU General Public
# License version 3. See COPYING for more information.

INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/include -I$(top_builddir)

pkgincludesubdir = $(includedir)/la-bencher
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = la-bencher.pc

SUBDIRS= . m4 src doc tests
OBJEXT=".lo .o"
ACLOCAL_AMFLAGS=-I m4

bin_SCRIPTS=la-bencher-config

EXTRA_DIST=PROJECTS

CLEANFILES=				\
	tags						\
	obj/*						\
	*~



distclean-local:
	rm -rf autom4te.cache
	rm -rf obj

DISTCLEANFILES= 												\
	_configs.sed													\
	config.sub														\
	config.guess													\
	configure.in													\
	configure															\
	ltmain.sh															\
	Makefile.in														\
	la-bencher-config la-bencher-config~	\
	include/la-bencher-config.h*					\
	include/stamp-h1a         						\
	doc/index.html												\
	obj/*																	\
	*~

ctags:
	ctags -R *

all-am: ctags
