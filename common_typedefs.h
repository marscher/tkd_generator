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

typedef std::vector<ug::vector3> CoordsArray;
typedef std::vector<int> IndexArray;

// subset indices
enum TKDSubsetTypes {
	LIPID = 0, CORNEOCYTE, BOUNDARY_CORN, BOUNDARY_LIPID
};

#endif /* COMMON_TYPEDEFS_H_ */
