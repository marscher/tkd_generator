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


typedef ug::vector3 v;
typedef std::vector<v> CoordsArray;
typedef std::vector<int> IndexArray;
typedef double number;


using ug::PI;
using std::endl;
using std::string;

#endif /* COMMON_TYPEDEFS_H_ */
