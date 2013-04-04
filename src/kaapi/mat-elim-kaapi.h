/**
 * \file   mat-elim-kaapi.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   April 2013
 * \brief  Header file for Gaussian Elimination using XKAAPI.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_KAAPI_H
#define F4RT_MAT_ELIM_KAAPI_H

#include <matrix.h>
#include "../mat-elim-tools.h"

void elimKAAPIC(Matrix& A, int blocksize);

void elimNaiveKAAPICModP1d(Matrix& A, int nthrds, int blocksize, uint64 prime);

void elimNaiveKAAPICModP1dPivot(Matrix& A, int nthrds, int blocksize, uint64 prime);
#endif
