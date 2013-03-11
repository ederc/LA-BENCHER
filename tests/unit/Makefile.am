# Copyright 2012 Christian Eder
# 
# This file is part of F4RT, licensed under the GNU General Public
# License version 3. See COPYING for more information.

INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/include -I$(top_srcdir)/src -I$(top_builddir)
AM_CPPFLAGS = -g3 -ggdb -Wall -D__F4RT_DEBUG -O0 -fomit-frame-pointer -pthread

BENCHMARK_CXXFLAGS = -O2 -pthread

SUBDIRS =

if ENABLE_DEBUG
LDADD = $(top_builddir)/src/.libs/libf4rt_debug.la
else
LDADD = $(top_builddir)/src/.libs/libf4rt.la
endif

# Put an entry for your test-mycomponent.C here (or in one of the other test groups).  
# Don't forget to add a test_mycomponent_SOURCES entry below.
#
# Once present, tests should remain present, being moved from one group to another.
# Note that among the available groups are NONCOMPILING_TESTS and FAILING_TESTS, which
# will be ignored by "make check" and can be temporary storage for problematic tests.
# a benchmarker, not to be included in check.
UNIT_TESTS =						\
	_add_mm128_epi32			\
	_mullo_mm128_epi32		\
	_dense_mat_transpose	\
	_dense_mat_copy				\
	_dense_mat_check

BENCHMARKS =            

EXTRA_PROGRAMS = $(NON_COMPILING_TESTS) $(BENCHMARKS)

check_PROGRAMS = $(UNIT_TESTS)

_add_mm128_epi32_SOURCES =			\
		_add_mm128_epi32.c
_mullo_mm128_epi32_SOURCES =		\
		_mullo_mm128_epi32.c
_dense_mat_transpose_SOURCES =	\
		_dense_mat_transpose.cpp
_dense_mat_copy_SOURCES =				\
		_dense_mat_copy.cpp
_dense_mat_check_SOURCES =			\
		_dense_mat_check.cpp

noinst_HEADERS =	
