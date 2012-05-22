/*
 * geometric_helper.h
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */

#ifndef GEOMETRIC_HELPER_H_
#define GEOMETRIC_HELPER_H_

#include "../common_typedefs.h"
#include "rotation_matrix.h"
namespace ug {
namespace tkd {

enum axis {
	xAxis = 0, yAxis, zAxis
};

enum {
	CLEAR=0
};

CoordsArray& operator<<(CoordsArray& array, const ug::vector3& vector);
CoordsArray& operator<<(CoordsArray& array, const int clear);

std::ostream& operator<<(std::ostream& out, const CoordsArray& arr);

ug::vector3 mirror(const ug::vector3& vec, const int axis);
CoordsArray mirror(const CoordsArray& coords, const int axis);

ug::vector3 translate(const ug::vector3& vec, const ug::vector3& offset);
CoordsArray translate(const CoordsArray&, const ug::vector3& offset);

//void reflect(ug::vector3& vOut, const ug::vector3& v, const ug::vector3& l, const number& c);

class myTransform {
public:
	myTransform(RotationMatrix& R, const ug::vector3& origin) :
		m_R(R), m_origin(origin) {};
	const ug::vector3 perform(const ug::vector3&) const;
	CoordsArray perform(const CoordsArray&);
protected:
	RotationMatrix& m_R;
	const ug::vector3& m_origin;
};

} // end of namespace tkd
} // end of namespace ug
#endif /* GEOMETRIC_HELPER_H_ */
