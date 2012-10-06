/*
 * vecComparator.h
 *
 *  Created on: 08.02.2012
 *      Author: marscher
 */

#ifndef VECCOMPARATOR_H_
#define VECCOMPARATOR_H_
#include "common/math/ugmath_types.h"
namespace ug {
namespace tkd {

struct vecComparator {

	// returns a < b
	bool operator()(const vector3&, const vector3&) const;
	
	static const number SMALL;
};

} /* namespace tkd */
} /* namespace ug */
#endif /* VECCOMPARATOR_H_ */
