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

//#include "matrix.h"

typedef unsigned long long uint64;
typedef unsigned int uint32;
typedef unsigned short uint16;
typedef unsigned char uint8;

typedef signed long long int64;
typedef signed int int32;
typedef signed short int16;
typedef signed char int8;

// matrix entry type:
// if uint16 entries are used we store uint64
// if float entries are used we store double
typedef uint64 mat; // mat entry type
typedef uint16 rmat; // real type of entry

#ifdef __cplusplus
extern "C" {
#endif
double countGEPFlops(uint32 m, uint32 n, uint64 prime);

//void cleanUpModP(Matrix& A, uint64 prime);

mat negInverseModP(mat a, uint64 p);
#ifdef __cplusplus
}
#endif
#endif
