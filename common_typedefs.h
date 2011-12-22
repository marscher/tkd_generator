/*
 * common_typedefs.h
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */

#ifndef COMMON_TYPEDEFS_H_
#define COMMON_TYPEDEFS_H_

#include <vector>
#include "common/math/ugmath_types.h"
#include "common/math/misc/math_util.h"

using ug::vector3;
using ug::matrix33;

using std::endl;
using std::string;

typedef std::vector<vector3> CoordsArray;
typedef std::vector<int> IndexArray;
typedef double number;

#endif /* COMMON_TYPEDEFS_H_ */
