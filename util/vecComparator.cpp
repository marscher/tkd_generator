/*
 * vecComparator.cpp
 *
 *  Created on: 30.03.2012
 *      Author: marscher
 */
#include "vecComparator.h"
namespace tkd {

// returns a < b
bool vecComparator::operator()(const ug::vector3& a,
		const ug::vector3& b) const {

	if (a.x < b.x - SMALL)
		return true;
	if (a.x > b.x + SMALL)
		return false;

	if (a.y < b.y - SMALL)
		return true;
	if (a.y > b.y + SMALL)
		return false;

	if (a.z < b.z - SMALL)
		return true;
	if (a.z > b.z + SMALL)
		return false;

	return false;
}

const number vecComparator::SMALL = 10E-6;

} /* namespace tkdGenerator */
