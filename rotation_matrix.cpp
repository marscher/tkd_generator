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
 * Rotate around y axis
 *
 +cos(t)  0  - sin(t)+
 |                   |
 |  0     1     0    |
 |                   |
 +sin(t)  0   cos(t) +
 */
void RotationMatrix::setAngle(const number& theta) {
	R[0][0] = cos(theta);
	R[0][2] = -sin(theta);
	R[2][0] = sin(theta);
	R[2][2] = cos(theta);
}

void RotationMatrix::setMirrorXAxis(const bool mirror) {
	if(mirror)
		R[1][1] = -1;
	else
		R[1][1] = 1;
}

RotationMatrix::RotationMatrix(const number& theta, const bool mirror) {
	R[0][0] = cos(theta);
	R[0][1] = 0;
	R[0][2] = -sin(theta);

	R[1][0] = 0;
	if(mirror)
		R[1][1] = -1;
	else
		R[1][1] = 1;
	R[1][2] = 0;

	R[2][0] = sin(theta);
	R[2][1] = 0;
	R[2][2] = cos(theta);
}

vector3 RotationMatrix::operator*(const vector3& v) {
	// temp values
	number x, y, z;
	x = R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2];
	y = R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2];
	z = R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2];

	return vector3(x, y, z);
}

} // end of namespace

