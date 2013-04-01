/**
 * \file   modular.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for modular tools.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#ifndef F4RT_MODULAR_H
#define F4RT_MODULAR_H

#include "matrix.h"

mat negInverseModP(mat a, uint64 p);
#endif
