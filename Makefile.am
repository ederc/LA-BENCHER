# Copyright 2012 Christian Eder
# 
# This file is part of F4RT, licensed under the GNU General Public
# License version 3. See COPYING for more information.

INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/include -I$(top_builddir)

pkgincludesubdir = $(includedir)/f4rt
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = f4rt.pc

SUBDIRS= . m4 src doc tests
OBJEXT=".lo .o"
ACLOCAL_AMFLAGS=-I m4

bin_SCRIPTS=f4rt-config

EXTRA_DIST=PROJECTS

CLEANFILES=				\
	tags						\
	*~

distclean-local:
	rm -rf autom4te.cache

DISTCLEANFILES= 									\
	_configs.sed										\
	config.sub											\
	config.guess										\
	configure.in										\
	configure												\
	ltmain.sh												\
	Makefile.in											\
	f4rt-config f4rt-config~				\
	include/f4rt-config.h*					\
	include/stamp-h1a         			\
	doc/index.html									\
	*~

ctags:
	ctags -R *

all-am: ctags
