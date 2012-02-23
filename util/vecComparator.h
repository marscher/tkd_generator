/*
 * vecComparator.h
 *
 *  Created on: 08.02.2012
 *      Author: marscher
 */

#ifndef VECCOMPARATOR_H_
#define VECCOMPARATOR_H_

namespace tkdGenerator {

struct vecComperator {

	// returns a < b
	bool operator()(const vector3& a, const vector3& b) const {
		
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
	
	static const number SMALL = 10E-6;
};


} /* namespace tkdGenerator */
#endif /* VECCOMPARATOR_H_ */
