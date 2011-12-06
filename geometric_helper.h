/*
 * geometric_helper.h
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */

#ifndef GEOMETRIC_HELPER_H_
#define GEOMETRIC_HELPER_H_

#include "common_typedefs.h"
#include "rotation_matrix.h"
namespace tkdGenerator {

void createPrism(vRef v1, vRef v2, vRef v3,
				 vRef v4, vRef v5, vRef v6,
				 CoordsArray& posOut, IndexArray& indsOut);

void createTetrahedron(vRef v1, vRef v2, vRef v3, vRef v4,
		CoordsArray& posOut, IndexArray& indsOut);

v myTransform(vRef v, RotationMatrix& R, vRef origin);

}

#endif /* GEOMETRIC_HELPER_H_ */
