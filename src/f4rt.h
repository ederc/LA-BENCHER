/**
 * \file   f4rt.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   February 2013
 * \brief  General header file for f4rt stuff.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "f4rt-config.h"

#include "matrix.h"
#include "mat-mult.h"
#include "mat-gen.h"

#define PACKAGE "F4RT"
#define VERSION "0.0.1"

void print_help(int exval);
