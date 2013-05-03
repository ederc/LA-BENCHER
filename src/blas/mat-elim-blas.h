/**
 * \file   mat-elim-seq.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for Gaussian Elimination using BLAS.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_BLAS_H
#define F4RT_MAT_ELIM_BLAS_H

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <f4rt-config.h>
#include "../mat-elim-tools.h"

#ifdef __F4RT_HAVE_LAPACK

#ifdef __cplusplus
extern "C" {
#endif
void elimBLAS(double *M, uint32 rows, uint32 cols, int nthrds, uint32 blocksize, uint64 prime);

#ifdef __cplusplus
}
#endif
#endif
#endif
