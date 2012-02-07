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
 * @param nOut stores normal vector of intersection line
 * @param p1 support vector of plane 1
 * @param n1 normal vector of plane 1
 * @param p2 support vector of plane 2
 * @param n2 normal vector of plane 1
 * @return true, if planes intersect; false if planes are parallel
 */
bool PlanePlaneIntersection3D(vector3& pOut, vector3& nOut, const vector3& p1,
		const vector3& n1, const vector3& p2, const vector3& n2) {

	// check planes are not parallel
	// which is the case iff | n1 * n2 | = 1
	if (fabs(VecDot(n1, n2)) - 1) < SMALL)
		return false;

	// intersection normal perpendicular to n1 and n2
	VecCross(nOut, n1, n2);

	// dertermine point on intersection line
	RayPlaneIntersection(pOut, p1, nOut, p2, n2);
	return true;
}

#endif /* PLANE_PLANE_INTERSECTION_H_ */
