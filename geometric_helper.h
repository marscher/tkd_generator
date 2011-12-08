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

//void createGeometricObject(int numVertices, const CoordsArray & posIn,
//		CoordsArray & posOut, IndexArray & indsOut);

CoordsArray& operator<<(CoordsArray& array, const v& vector);

std::ostream& operator<<(std::ostream& out, const CoordsArray& arr);

class myTransform {
public:
	myTransform(RotationMatrix&, const v& origin);
	const v perform(const v&);
	CoordsArray perform(CoordsArray&);
protected:
	RotationMatrix& m_R;
	const v& m_origin;
};

}

#endif /* GEOMETRIC_HELPER_H_ */
