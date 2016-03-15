/*
 * vecComparator.cpp
 *
 *  Created on: 30.03.2012
 *      Author: marscher
 */
#include "./vecComparator.h"
namespace std {
/**
 * To compare vector3 with < operator (lexicographically sorted by component).
 * Use with care!
 */
bool less<vector3>::operator()(const vector3& a, const vector3& b) const {
	return a < b;
}

/**
 * assumes a == b iff !(a < b) && !(b < a)
 */
bool equal_to<vector3>::operator()(const vector3& a, const vector3& b) const {
	return a == b;
}

} // end of namespace std

// ug exposed comparator
namespace ug {

bool operator==(const vector3& a, const vector3& b) {
	return !(a < b) && !(b < a);
}

} /* namespace ug */
