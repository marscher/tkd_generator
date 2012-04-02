/*
 * geometric_helper.cpp
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */
#include "geometric_helper.h"
#include "common/log.h"
#include "common/math/math_vector_matrix/math_matrix.h"
#include "common/math/math_vector_matrix/math_matrix_vector_functions.h"
namespace tkd {

/**
 * This method first rotates vector v with rotation matrix R and then translates
 * the result by origin:
 * v' = (R*v)+origin
 * note: v' is a copy. The transformation creates a new instance
 * @param v vector to transform
 */
const ug::vector3 myTransform::perform(const ug::vector3& v) const {
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
CoordsArray& operator<<(CoordsArray& array, const ug::vector3& vector) {
	array.push_back(vector);
	return array;
}

/**
 * clears the vector in operator << style
 */
CoordsArray& operator<<(CoordsArray& array, const int clear) {
	array.clear();
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

ug::vector3 mirror(const ug::vector3& vec, const int axis) {
	ug::vector3 result;
	result.x = vec.x * (axis == xAxis ? -1 : 1);
	result.y = vec.y * (axis == yAxis ? -1 : 1);
	result.z = vec.z * (axis == zAxis ? -1 : 1);
	return result;
}

CoordsArray mirror(const CoordsArray& coords, const int axis) {
	CoordsArray result;
	for (size_t i = 0; i < coords.size(); i++) {
		result.push_back(mirror(coords[i], axis));
	}
	return result;
}

ug::vector3 translate(const ug::vector3& vec, const ug::vector3& offset) {
	ug::vector3 result;
	result.x = vec.x + offset.x;
	result.y = vec.y + offset.y;
	result.z = vec.z + offset.z;
	return result;
}

CoordsArray translate(const CoordsArray& coords, const ug::vector3& offset) {
	CoordsArray result;
	for (size_t i = 0; i < coords.size(); i++) {
		result.push_back(translate(coords[i], offset));
	}
	return result;
}

//void reflect(ug::vector3& vOut, const ug::vector3& v, const ug::vector3& a, const number& c) {
//	number s = 2 * (VecDot(v, a) - c) / VecDot(a, a);
//	VecScale(vOut, a, s);
//	VecSubtract(vOut, v, vOut);
//}

} // end of namespace tkdGenerator
