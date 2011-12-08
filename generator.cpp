/*
 * generator.cpp
 *
 *  Created on: 08.12.2011
 *      Author: marscher
 */
#include "generator.h"
namespace tkdGenerator {

void Generator::createTop(const v& offset) {

	myTransform t(R, offset);
	CoordsArray obenInnen;
	obenInnen << v1 << v2 << v3 << v4 << v5 << v6;

	// G(Ki -> ObenAussenPr)
	CoordsArray obenAussenPr;
	obenAussenPr << v3 << v6 << v7 << v8 << v9 << v10;

	CoordsArray obenAussenPr_leftT;
	obenAussenPr_leftT << v3 << v7 << v11 << v6;

	CoordsArray obenAussenPr2T_prism;
	obenAussenPr2T_prism << v3 << v6 << v11 << v2 << v5 << v12;

	CoordsArray obenAussenPr_righT;
	obenAussenPr_righT << v2 << v12 << v13 << v5;

	// create G(Ki -> ObenInnen) = 6 prism with equilateral sites
	// step angle 60° = PI/3 (rad)
	for (number r = 0; r < 2 * PI; r += PI / 3) {
		// sets the angle of rotation matrix which is applied on all points in this loop
		R.setAngle(r);
		CoordsArray c = t.perform(obenInnen);
		createGeometricObject(c);
	}

	// create G(Ki -> ObenAussenPr2T) and G(Ki -> ObenAussenPr
	for (number r = 2. / 3 * PI; r < 2 * PI; r += (2. / 3 * PI)) {
		// sets the angle of rotation matrix which is applied on all points in this loop
		R.setAngle(r);

		// create G(Ki -> ObenAussenPr):
		createGeometricObject(t.perform(obenAussenPr));

		/* tetrahedron left of prism of ObenAussenPr2T which
		 shares 3 points with prism of ObenAussenPr */
		createGeometricObject(t.perform(obenAussenPr_leftT));

		// Prism of ObenAussenPr2T
		createGeometricObject(t.perform(obenAussenPr2T_prism));

		/* tetrahedron right of prism ObenAussenPr2T which shares
		 3 points with prism obenAussenPr2T_prism */
		createGeometricObject(t.perform(obenAussenPr_righT));
	}
}

void Generator::createMiddle(const v& origin) {
	myTransform t(R, origin);
	// create G(Ki -> ObenInnen) = 6 prism with equilateral sites
	// step angle 60° = PI/3 (rad)

	CoordsArray obenInnen;
	obenInnen << v1 << v2 << v3 << v4 << v5 << v6;
	for (number r = 0; r < 2 * PI; r += PI / 3) {
		// sets the angle of rotation matrix which is applied on all points in this loop
		R.setAngle(r);
		createGeometricObject(t.perform(obenInnen));
	}
}

/**
 * creates a tkd starting construction in (0, 0, 0)
 */
void Generator::createTKD() {
	createTKD(origin);
}

/**
 * creates a tkd starting with
 */
void Generator::createTKD(const v& offset) {
	createTop(offset);
	//TODO does this work if offset is negative?
	v offset_h(offset.x, offset.y + (-h), offset.z);

	createMiddle(offset_h);

	// offset is same like middle because bottom is mirrored around x axis
	R.setMirrorXAxis(true);
	// in fact this is the bottom part
	createTop(offset_h);
}

/**
 * init basis geometrical parameters
 */
void Generator::init() {
	// global index for grid vertices
	index = 0;

	s = 1 / sqrt(3) * (w - 2 * a);
	// height of base triangle of top inner prism
	g = sqrt(3) * a / 2;
	// height of base triangle of tetrahedron of ObenAussenPr2T
	b = sqrt(3) * s / 2;
	// height of hyphotenuse s in base triangle of tetrahedron of ObenAussenPr2T
	// describing offset in x direction
	hs = sqrt(3) * s / 4;
	// segment of hypothenuse s describing offset in x direction
	q = s / 4;

	UG_LOG("g: " << g << endl);
	UG_LOG("s: " << s << endl);
	UG_LOG("b: " << b << endl);

	// assign common vertices needed for construction.
	// tODO can this somehow be const?
	v1 = v(0, 0, 0);
	v2 = v(a, 0, 0);
	v3 = v(a / 2, 0, g);
	v4 = v(0, h, 0);
	v5 = v(a, h, 0);
	v6 = v(a / 2, h, g);
	v7 = v(a / 2, 0, g + s);
	v8 = v(-a / 2, 0, g);
	v9 = v(-a / 2, h, g);
	v10 = v(-a / 2, 0, g + s);
	v11 = v(a / 2 + hs, 0, g + q);
	v12 = v(a + hs, 0, q);
	v13 = v(a + hs + q, 0, -q);
}

void Generator::createGeometricObject(const CoordsArray & posIn) {
	//TODO determine which call has empty posIn!
	if(posIn.empty())
		return;

	uint numCoords = posIn.size();

	switch (numCoords) {
	case Tetrahedron:
		break;
	case Pyramid:
		break;
	case Prism:
		break;
	case Hexahedron:
		break;
	default:
		UG_LOG(" numcoords ist kein geometrischer typ: " << numCoords);
		assert(false);
	}

	// insert contents of posIn at end of global position array
	posOut.insert(posOut.end(), posIn.begin(), posIn.end());

	// tell lib_grid how much vertex indices are belonging to this geometric object
	indsOut.push_back(posIn.size());
	// enum each node of this geometric object to assign unique node index
	for (size_t current = index; index < current + numCoords; index++) {
		indsOut.push_back(index);
	}
}

} // end of namespace tkdGenerator
