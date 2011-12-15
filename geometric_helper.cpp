/*
 * geometric_helper.cpp
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */
#include "geometric_helper.h"
#include "common/log.h"

namespace tkdGenerator {

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

CoordsArray myTransform::perform(const CoordsArray& coords) {
	CoordsArray result;
	for (size_t i = 0; i < coords.size(); i++) {
		result.push_back(perform(coords[i]));
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
std::ostream & operator <<(std::ostream & out, const CoordsArray& coords) {
	for (size_t i = 0; i < coords.size(); i++) {
		out << coords[i] << ", ";
	}
	return out;
}

v mirror(const v& vec, const int axis) {
	v result;
	result.x = vec.x * (axis == xAxis ? -1 : 1);
	result.y = vec.y * (axis == yAxis ? -1 : 1);
	result.z = vec.z * (axis == zAxis ? -1 : 1);
	return result;
}

CoordsArray mirror(const CoordsArray coords, const int axis) {
	CoordsArray result;
	for (size_t i = 0; i < coords.size(); i++) {
		result.push_back(mirror(coords[i], axis));
	}
	return result;
}

v translate(const v& vec, const v& offset) {
	v result;
	result.x = vec.x + offset.x;
	result.y = vec.y + offset.y;
	result.z = vec.z + offset.z;
	return result;
}

CoordsArray translate(const CoordsArray& coords, const v& offset) {
	CoordsArray result;
	for (size_t i = 0; i < coords.size(); i++) {
		result.push_back(translate(coords[i], offset));
	}
	return result;
}

} // end of namespace tkdGenerator
