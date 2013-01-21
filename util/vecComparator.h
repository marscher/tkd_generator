/*
 * vecComparator.h
 *
 *  Created on: 08.02.2012
 *      Author: marscher
 */

#ifndef VECCOMPARATOR_H_
#define VECCOMPARATOR_H_

#include "common/math/ugmath_types.h"
#include "common/math/misc/math_util.h"

namespace ug {

template <std::size_t dim, typename T>
bool operator<(const MathVector<dim, T>& a, const MathVector<dim, T>& b) {
	for (int i = 0; i < dim; i++) {
		if(a.m_data[i] < b.m_data[i] - SMALL)
			return true;
		if(a.m_data[i] > b.m_data[i] + SMALL)
			return false;
	}

	return false;
}

// already template defined for operator==, so specialize it only for 3d
bool operator==(const vector3& a, const vector3& b) {
	return !(a < b) && !(b < a);
}

} // end of namespace ug

namespace std {
using ug::vector3;

/**
 * To compare vector3 with < operator (lexicographically sorted by component).
 * Use with care!
 */
template<> struct less<vector3> : public binary_function<vector3, vector3, bool> {
	bool operator()(const vector3& a, const vector3& b) const {
		return a < b;
	}
};

/**
 * assumes a == b iff !(a < b) && !(b < a)
 */
template<> struct equal_to<vector3> : public binary_function<vector3, vector3,
		bool> {
	bool operator()(const vector3& a, const vector3& b) const {
		return a == b;
	}
};

} // end of namespace std

#endif /* VECCOMPARATOR_H_ */
