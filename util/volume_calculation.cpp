/*
 * volume_calculation_impl.h
 *
 *  Created on: 11.01.2012
 *      Author: marscher
 */

#ifndef VOLUME_CALCULATION_IMPL_H_
#define VOLUME_CALCULATION_IMPL_H_

#include "volume_calculation.h"
//#include "lib_grid/lg_base.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"
#include "lib_algebra/common/operations_vec.h"

#include <vector>
#include <algorithm>
#include <iterator>
using namespace std;

namespace ug {

number CalculateVolume(const Volume& vol) {
	number result = -1;
	switch (vol.reference_object_id()) {
	case ROID_TETRAHEDRON:
		result = CalculateVolume(static_cast<Tetrahedron>(vol));
		break;
	case ROID_PRISM:
		result = CalculateVolume(static_cast<Prism>(vol));
		break;
	case ROID_PYRAMID:
		result = CalculateVolume(static_cast<Pyramid>(vol));
		break;
	case ROID_HEXAHEDRON:
		result = CalculateVolume(static_cast<Hexahedron>(vol));
		break;
	default:
		break;
	}

	return result;
}

number CalculateVolume(const Tetrahedron& tet) {
	number result;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	vector3& a = aaPos[tet.vertex(0)];
	vector3& b = aaPos[tet.vertex(1)];
	vector3& c = aaPos[tet.vertex(2)];
	vector3& d = aaPos[tet.vertex(3)];
	vector3 ad, bd, cd;

	VecSubtract(ad, a, d);
	VecSubtract(bd, b, d);
	VecSubtract(cd, c, d);

	vector3 cross;
	VecCross(cross, bd, cd);
	result = VecProd(cross, ad) / 6;
	return result;
}

number CalculateVolume(const Prism& prism) {
	Grid grid;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	const VolumeVertices* cvolvert = &prism;
	VolumeVertices* volvert = const_cast<VolumeVertices*>(cvolvert);
	// fixme compiler tries to cast volvert to FaceVertices*
	// (perhaps because overrided template CalculateCenter() is instanced in this unit another time with FaceVertices*
	vector3 centerPos = CalculateCenter(volvert, aaPos);

	VertexBase* center = *grid.create<Vertex>();

	aaPos[center] = centerPos;

	VertexBase* v0 = prism.vertex(0);
	VertexBase* v1 = prism.vertex(1);
	VertexBase* v2 = prism.vertex(2);
	VertexBase* v3 = prism.vertex(3);
	VertexBase* v4 = prism.vertex(4);
	VertexBase* v5 = prism.vertex(5);

	// top
	TetrahedronDescriptor t1(v0, v1, v2, center);

	PyramidDescriptor p1(v0, v1, v4, v3, center);
	PyramidDescriptor p2(v1, v2, v5, v4, center);
	PyramidDescriptor p3(v0, v3, v5, v2, center);

	// bottom
	TetrahedronDescriptor t2(v3, v4, v5, center);

	grid.create<Pyramid>(p1);
	grid.create<Pyramid>(p2);
	grid.create<Pyramid>(p3);

	grid.create<Tetrahedron>(t1);
	grid.create<Tetrahedron>(t2);

	return CalculateVolume(grid.begin<Volume>(), grid.end<Volume>());
}

number CalculateVolume(const Pyramid& pyramid) {
	number result;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	vector3& a = aaPos[pyramid.vertex(0)];
	vector3& b = aaPos[pyramid.vertex(1)];
	vector3& c = aaPos[pyramid.vertex(2)];
	vector3& d = aaPos[pyramid.vertex(3)];
	vector3& top = aaPos[pyramid.vertex(4)];

	//TODO is face[0] base area?
	FaceDescriptor base = pyramid.face(0);
	vector3 center = CalculateCenter(dynamic_cast<FaceVertices*>(&base), aaPos);
	number h = VecLength(top -= center);
	vector3 cross1, cross2;
	VecCross(cross1, a, b);
	VecCross(cross2, c, d);
	number A = VecLength(cross1) - VecLength(cross2);
	result = 1.0 / 3 * A * h;

	return result;
}

/**
 * Algorithm (14) from J.Grandy "Efficient Computation of Volume of Hexahedral Cells"
 */
number CalculateVolume(const Hexahedron& hexa) {
	number result;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	// bottom quad
	vector3& a = aaPos[hexa.vertex(0)];
	vector3& b = aaPos[hexa.vertex(1)];
	vector3& c = aaPos[hexa.vertex(2)];
	vector3& d = aaPos[hexa.vertex(3)];
	// top quad
	vector3& e = aaPos[hexa.vertex(4)];
	vector3& f = aaPos[hexa.vertex(5)];
	vector3& g = aaPos[hexa.vertex(6)];
	vector3& h = aaPos[hexa.vertex(7)];

	// determine long diagonal
	vector<number> diagonalLength;
	vector<number>::iterator iter;

	// diagonals
	vector3 ag, bh, ce, df;
	// longest diagonal
	vector3* LD = NULL;
	VecSubtract(ag, a, g);
	VecSubtract(bh, b, h);
	VecSubtract(ce, c, e);
	VecSubtract(df, d, f);

	diagonalLength.push_back(VecLength(ag));
	diagonalLength.push_back(VecLength(bh));
	diagonalLength.push_back(VecLength(ce));
	diagonalLength.push_back(VecLength(df));
	// determine longest diagonal
	iter = max_element(diagonalLength.begin(), diagonalLength.end());
	uint pos = distance(diagonalLength.begin(), iter);
	switch (pos) {
	case 0:
		LD = &ag;
		break;
	case 1:
		LD = &bh;
		break;
	case 2:
		LD = &ce;
		break;
	case 3:
		LD = &df;
	}

	// matrices to calculate determinant
	matrix33 d1, d2, d3;
	vector3 LDa, ba, ea, fg, ca, gd;

	VecSubtract(LDa, *LD, a);
	VecSubtract(ba, b, a);
	VecSubtract(ea, e, a);
	VecSubtract(fg, f, g);
	VecSubtract(ca, c, a);
	VecSubtract(gd, g, d);

	d1.assign(LDa, 0);
	d1.assign(ba, 1);
	d1.assign(df, 2);

	d2.assign(LDa, 0);
	d2.assign(ea, 1);
	d2.assign(fg, 2);

	d3.assign(LDa, 0);
	d3.assign(ca, 1);
	d3.assign(gd, 2);

	return result = Determinant(d1) + Determinant(d2) + Determinant(d3);
}

number CalculateVolume(geometry_traits<Volume>::iterator begin, geometry_traits<Volume>::iterator end) {
	number result = 0;
	for (VolumeIterator iter = begin; iter != end; iter++) {
		result += CalculateVolume(**iter);
	}
	return result;
}

} // end of namespace ug
#endif /* VOLUME_CALCULATION_IMPL_H_ */
