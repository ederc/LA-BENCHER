/**
 * \file   mat-elim-omp.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for Gaussian Elimination using OMP.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_OMP_H
#define F4RT_MAT_ELIM_OMP_H

#include <matrix.h>
#include "../mat-elim-tools.h"

void elimOMP(Matrix& A, int blocksize);
void elimNaiveOMPModP1dOuter(Matrix& A, int nthrds, int blocksize, uint64 prime);
#endif
