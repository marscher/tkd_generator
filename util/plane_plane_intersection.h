/*
 * plane_plane_intersection.h
 *
 *  Created on: 07.02.2012
 *      Author: marscher
 */

#ifndef PLANE_PLANE_INTERSECTION_H_
#define PLANE_PLANE_INTERSECTION_H_
#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "common/math/misc/math_util.h"

/**
 * @param pOut stores support vector of intersection line
 * @param dOut stores the direction vector of the intersection line
 * @param p1 support vector of plane 1
 * @param n1 normal vector of plane 1
 * @param p2 support vector of plane 2
 * @param n2 normal vector of plane 1
 * @return true, if planes intersect; false if planes are parallel
 */
bool PlanePlaneIntersection(vector3& pOut, vector3& dOut, const vector3& p1,
		const vector3& n1, const vector3& p2, const vector3& n2) {

	// check planes are not parallel
	// which is the case iff | n1 * n2 | = 1
	vector3 tn1, tn2;

	VecNormalize(tn1, n1);
	VecNormalize(tn2, n2);

	if (fabs(VecDot(tn1, tn2)) - 1 < SMALL)
		return false;

	// intersection normal perpendicular to n1 and n2
	VecCross(dOut, n1, n2);

	//	use this temporary direction to find pOut.
	vector3 tdir;
	VecCross(tdir, dOut, n1);

	// dertermine point on intersection line
	number tmp;
	return RayPlaneIntersection(pOut, tmp, p1, tdir, p2, n2);
}

#endif /* PLANE_PLANE_INTERSECTION_H_ */
