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

typedef ug::vector3 v;
typedef std::vector<v> CoordsArray;
typedef std::vector<int> IndexArray;
typedef const v& vRef;
typedef double number;

using std::endl;
using std::string;

#endif /* COMMON_TYPEDEFS_H_ */
