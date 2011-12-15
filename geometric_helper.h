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

enum axis {
	xAxis = 0, yAxis = 1, zAxis = 2
};

CoordsArray& operator<<(CoordsArray& array, const v& vector);
std::ostream& operator<<(std::ostream& out, const CoordsArray& arr);

v mirror(const v& vec, const int axis);
CoordsArray mirror(const CoordsArray coords, const int axis);

v translate(const v& vec, const v& offset);
CoordsArray translate(const CoordsArray&, const v& offset);

class myTransform {
public:
	myTransform(RotationMatrix& R, const v& origin) :
		m_R(R), m_origin(origin) {};
	const v perform(const v&);
	CoordsArray perform(const CoordsArray&);
protected:
	RotationMatrix& m_R;
	const v& m_origin;
};

}

#endif /* GEOMETRIC_HELPER_H_ */
