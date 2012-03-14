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

namespace tkd {

enum axis {
	xAxis = 0, yAxis, zAxis
};

CoordsArray& operator<<(CoordsArray& array, const vector3& vector);
std::ostream& operator<<(std::ostream& out, const CoordsArray& arr);

vector3 mirror(const vector3& vec, const int axis);
CoordsArray mirror(const CoordsArray& coords, const int axis);

vector3 translate(const vector3& vec, const vector3& offset);
CoordsArray translate(const CoordsArray&, const vector3& offset);

//void reflect(vector3& vOut, const vector3& v, const vector3& l, const number& c);

class myTransform {
public:
	myTransform(RotationMatrix& R, const vector3& origin) :
		m_R(R), m_origin(origin) {};
	const vector3 perform(const vector3&);
	CoordsArray perform(const CoordsArray&);
protected:
	RotationMatrix& m_R;
	const vector3& m_origin;
};

}

#endif /* GEOMETRIC_HELPER_H_ */
