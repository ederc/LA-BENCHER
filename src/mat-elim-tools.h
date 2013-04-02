/**
 * \file   mat-elim-tools.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for GEP tools.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_TOOLS_H
#define F4RT_MAT_ELIM_TOOLS_H

#include "matrix.h"

double countGEPFlops(uint32 m, uint32 n);

void cleanUpModP(Matrix& A, uint64 prime);

mat negInverseModP(mat a, uint64 p);
#endif
