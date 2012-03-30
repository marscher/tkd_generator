/*
 * rotation_matrix.cpp
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */

#include "rotation_matrix.h"
#include "common/math/math_vector_matrix/math_matrix_vector_functions.h"
// needed for M_PI macro
#define _USE_MATH_DEFINES
#include <cmath>

namespace tkd {

/**
 * sets rotation angle around z axis in degree
 */
void RotationMatrix::setAngle(const number& deg) {
	number theta = deg * M_PI / 180.0;

	m_data[0][0] = cos(theta);
	m_data[0][1] = -sin(theta);
	m_data[1][0] = sin(theta);
	m_data[1][1] = cos(theta);
}

/**
 * sets whether to rotate around z axis
 * @param mirror mirror or not?
 */
void RotationMatrix::setMirrorZAxis(const bool mirror) {
	if (mirror)
		m_data[2][2] = -1;
	else
		m_data[2][2] = 1;
}

/**
 * @param deg initial rotatation angle in degree
 * @param mirror true if mirroring around z axis should be applied
 */
RotationMatrix::RotationMatrix(const number& deg, const bool mirror) {
	number theta = deg * M_PI / 180.0;

	m_data[0][0] = cos(theta);
	m_data[0][1] = -sin(theta);
	m_data[0][2] = 0;

	m_data[1][0] = sin(theta);
	m_data[1][1] = cos(theta);
	m_data[1][2] = 0;

	m_data[2][0] = 0;
	m_data[2][1] = 0;

	if (mirror)
		m_data[2][2] = -1;
	else
		m_data[2][2] = 1;
}

/**
 * performs rotation and, if set, mirroring around z axis.
 * @param vec vector to rotate (and mirror)
 */
vector3 RotationMatrix::operator*(const vector3& vec) const {
	vector3 result;
	MatVecMult(result, *this, vec);
	return result;
}

} // end of namespace
