/*
 * common_typedefs.h
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */

#ifndef COMMON_TYPEDEFS_H_
#define COMMON_TYPEDEFS_H_

#include <vector>
#include <common/types.h>
#include "common/math/ugmath_types.h"
#include "common/math/misc/math_util.h"

using ug::vector3;
using ug::matrix33;

using std::endl;
using std::string;
using std::vector;

typedef vector<vector3> CoordsArray;
typedef vector<int> IndexArray;

// subset indices
enum Subset {
	LIPID = 0, CORNEOCYTE
};

#endif /* COMMON_TYPEDEFS_H_ */
