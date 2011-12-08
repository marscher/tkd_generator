/*
 * geometric_helper.cpp
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */
#include "geometric_helper.h"
#include <cassert>
#include "common/log.h"
#include "common/assert.h"

namespace tkdGenerator {

myTransform::myTransform(RotationMatrix& R, const v& origin) :
		m_R(R), m_origin(origin) {
}

/**
 * This method first rotates vector v with rotation matrix R and then translates
 * the result by origin:
 * v' = (R*v)+origin
 * note: v' is a copy. The transformation creates a new instance
 * @param v vector to transform
 */
const v myTransform::perform(const v& v) {
	return (m_R * v) += m_origin;
}

CoordsArray myTransform::perform(CoordsArray& coords) {
	CoordsArray result(coords);
	for (size_t i = 0; i < coords.size(); i++) {
		result[i] = perform(coords[i]);
	}
	return result;
}

/**
 * simply pushes given vector at back of given CoordsArray
 */
CoordsArray& operator<<(CoordsArray& array, const v& vector) {
	array.push_back(vector);
	return array;
}

/**
 * print vector of vectors
 */
std::ostream& operator<<(std::ostream& out, std::vector<v>& arr) {
	for (size_t i = 0; i < arr.size(); i++) {
		out << arr[i] << endl;
	}
	return out;
}

} // end of namespace tkdGenerator
