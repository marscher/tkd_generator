/*
 * rotation_matrix.cpp
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */

#include "rotation_matrix.h"
#include <cmath>

namespace tkdGenerator {

/**
 * sets rotation angle around z axis in degree
 */
void RotationMatrix::setAngle(const number& deg) {
	number theta = deg * PI / 180.0;

	R[0][0] = cos(theta);
	R[0][1] = -sin(theta);
	R[1][0] = sin(theta);
	R[1][1] = cos(theta);
}

/**
 * sets whether to rotate around z axis
 * @param mirror mirror or not?
 */
void RotationMatrix::setMirrorZAxis(const bool mirror) {
	if (mirror)
		R[2][2] = -1;
	else
		R[2][2] = 1;
}

/**
 * @param deg initial rotatation angle in degree
 * @param mirror true if mirroring around z axis should be applied
 */
RotationMatrix::RotationMatrix(const number& deg, const bool mirror) {
	number theta = deg * PI / 180.0;

	R[0][0] = cos(theta);
	R[0][1] = -sin(theta);
	R[0][2] = 0;

	R[1][0] = sin(theta);
	R[1][1] = cos(theta);
	R[1][2] = 0;

	R[2][0] = 0;
	R[2][1] = 0;

	if (mirror)
		R[2][2] = -1;
	else
		R[2][2] = 1;
}

/**
 * performs rotation and if set mirroring around z axis. Result is a copy
 * @param vec vector to rotate (and mirror)
 */
v RotationMatrix::operator*(const v& vec) {
	v result;

	result.x = R[0][0] * vec[0] + R[0][1] * vec[1] + R[0][2] * vec[2];
	result.y = R[1][0] * vec[0] + R[1][1] * vec[1] + R[1][2] * vec[2];
	result.z= R[2][0] * vec[0] + R[2][1] * vec[1] + R[2][2] * vec[2];

	return result;
}

} // end of namespace
