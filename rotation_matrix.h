/*
 * rotation_matrix.h
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */

#include "common_typedefs.h"

#ifndef _ROTATION_MATRIX_
#define _ROTATION_MATRIX_

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
class RotationMatrix {
public:
	void setMirrorZAxis(const bool);
	void setAngle(const number& theta);
	RotationMatrix(const number& theta = 0, const bool mirror = false);
	v operator*(const v& vector);

protected:
	number R[3][3];
};

} // end of namespace
#endif
