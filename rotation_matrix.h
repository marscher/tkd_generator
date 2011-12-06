/*
 * rotation_matrix.h
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */

#include "common_typedefs.h"

#ifndef _ROTATION_MATRIX_
#define _ROTATION_MATRIX_

using ug::vector3;

namespace tkdGenerator {

/*
 +cos(t)  0  - sin(t)+
 |                   |
 |  0     1     0    |
 |                   |
 +sin(t)  0   cos(t) +
 */
class RotationMatrix {
public:
	void setMirrorXAxis(const bool);
	void setAngle(const number& theta);
	RotationMatrix(const number& theta, const bool mirror = false);
	vector3 operator*(const vector3& v);

protected:
	number R[3][3];
};

} // end of namespace
#endif
