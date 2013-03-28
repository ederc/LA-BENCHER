/**
 * \file   mat-elim-seq.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for sequential Gaussian Elimination.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MAT_ELIM_SEQ_H
#define F4RT_MAT_ELIM_SEQ_H

#include <matrix.h>

void elimSEQ(Matrix& A, int blocksize);
#endif
