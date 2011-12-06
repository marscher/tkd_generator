/*
 * geometric_helper.cpp
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */
#include "geometric_helper.h"

namespace tkdGenerator {

// global index for nodes
static unsigned int index = 0;

void createPrism(vRef v1, vRef v2, vRef v3, vRef v4, vRef v5,
		vRef v6, CoordsArray& posOut, IndexArray& indsOut) {
	posOut.push_back(v1);
	posOut.push_back(v2);
	posOut.push_back(v3);

	posOut.push_back(v4);
	posOut.push_back(v5);
	posOut.push_back(v6);

	// 6 nodes
	indsOut.push_back(6);
	// enum each node of prism to assign unique node index
	for (unsigned int current = index; index < current + 6; index++) {
		indsOut.push_back(index);
	}
}

/**
 * please be sure to pass the vertices in the correct order:
 * v1, v2, v3: bottom-vertices in counterclockwise order (if viewed from the top).
 * v4: top
 */
void createTetrahedron(vRef v1, vRef v2, vRef v3, vRef v4,
		CoordsArray& posOut, IndexArray& indsOut) {
	posOut.push_back(v1);
	posOut.push_back(v2);
	posOut.push_back(v3);
	// top
	posOut.push_back(v4);

	// 4 nodes
	indsOut.push_back(4);
	// enum each node of prism to assign unique node index
	for (unsigned int current = index; index < current + 4; index++) {
		indsOut.push_back(index);
	}
}

/**
 * This method first rotates vector v with rotation matrix R and then translates
 * the result by origin:
 * v' = (R*v)+origin
 * @param v vector to translate
 * @param R rotation matrix to use
 * @param origin offset with which Rv is translated
 */
v myTransform(vRef v, RotationMatrix& R, vRef origin) {
	return (R * v) += origin;
}

} // end of namespace tkdGenerator
