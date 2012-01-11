/*
 * rotation_matrix.h
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */

#ifndef _ROTATION_MATRIX_
#define _ROTATION_MATRIX_

#include "common_typedefs.h"

namespace tkdGenerator {

/**
 * class to perform transformation of rotation around z axis with following matrix
 *
 +cos(r)  - sin(r)  0+
 |                   |
 |sin(r)   cos(r)   0|
 |                   |
 +  0        0      1+
 */
class RotationMatrix : public matrix33 {
public:
	void setMirrorZAxis(const bool);
	void setAngle(const number& theta);
	RotationMatrix(const number& theta = 0, const bool mirror = false);
	vector3 operator*(const vector3& vector) const;
};

} // end of namespace
#endif
