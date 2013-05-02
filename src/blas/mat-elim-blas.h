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

#ifdef __F4RT_HAVE_OPENBLAS
#include <matrix.h>
#include "../mat-elim-tools.h"

// HACK FOR "template with C linkage error" from LAPACK
#include <complex.h>
#define lapack_complex_float float _Complex
#define lapack_complex_double double _Complex

#include <lapacke.h> // from OpenBLAS

#ifdef __cplusplus
extern "C" {
#endif
void elimBLAS(Matrix& A, uint32 blocksize, uint64 prime);

#ifdef __cplusplus
}
#endif
#endif
#endif
