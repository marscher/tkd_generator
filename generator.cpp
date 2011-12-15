/*
 * generator.cpp
 *
 *  Created on: 08.12.2011
 *      Author: marscher
 */
#include "generator.h"
namespace tkdGenerator {

/**
 * creates top and bottom part of tkd with given offset and rotationOffset
 * @param offset
 * @param rotationOffset
 */
void Generator::createTop(const v& offset, const number rotationOffset) {
	myTransform t(R, offset);

	// create G(Ki -> ObenInnen) = 6 prism with equilateral sites
	// step angle 60°
	for (int r = rotationOffset; r < 360 + rotationOffset; r += 60) {
		// sets the angle of rotation matrix which is applied on all points in this loop
		R.setAngle(r);
		createGeometricObject(t.perform(obenInnen));
	}

	//  create G(Ki -> ObenAussenPr)
	for (int r = rotationOffset; r < 360; r += 120) {
		// sets the angle of rotation matrix which is applied on all points in this loop
		R.setAngle(r);

		// create G(Ki -> ObenAussenPr):
		createGeometricObject(t.perform(obenAussenPrism));
	}

	// create G(Ki -> ObenAussenPr2T) with offset of 60° to ObenAussenPr
	for (int r = rotationOffset + 60; r < 360 + rotationOffset; r += 120) {
		// Prism of ObenAussenPr2T
		R.setAngle(r);
		createGeometricObject(t.perform(obenAussenPr2T_prism));

		/* tetrahedron left of prism of ObenAussenPr2T which
		 shares 3 points with prism of ObenAussenPr */
		createGeometricObject(t.perform(obenAussenPr_leftTetrahedrson));

		/* tetrahedron right of prism ObenAussenPr2T which shares
		 3 points with prism obenAussenPr2T_prism */
		createGeometricObject(t.perform(obenAussenPr_rightTetrahedron));
	}
}

void Generator::createMiddle(const v& origin) {
	myTransform t(R, origin);

	for (int r = 0; r < 360; r += 60) {
		// sets the angle of rotation matrix which is applied on all points in this loop
		R.setAngle(r);
		createGeometricObject(t.perform(obenInnen));
	}

	for (int r = 0; r < 360; r += 120) {
		R.setAngle(r);
		createGeometricObject(t.perform(mitteAussenHexahedron));
		createGeometricObject(t.perform(mitteAussenH2Pr_pyramid));
		createGeometricObject(t.perform(mitteAussen2PrH_pyramid));
	}

	for (int r = 60; r < 360; r += 120) {
		R.setAngle(r);

		createGeometricObject(t.perform(mitteAussenP1));
		createGeometricObject(t.perform(mitteAussenP2));

		createGeometricObject(t.perform(mitteAussenH2Pr_tetrahedron));
		createGeometricObject(t.perform(mitteAussen2PrH_tetrahedron));
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
	createTop(offset, 0);
	v offset_h(offset.x, offset.y, offset.z - h);

	createMiddle(offset_h);

	// offset is same like middle because bottom is mirrored around x axis
	R.setMirrorZAxis(true);
	// in fact this is the bottom part, rotated with 60° relative to top
	createTop(offset_h, 60);
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

	UG_LOG("g: " << g << endl);
	UG_LOG("s: " << s << endl);
	UG_LOG("b: " << b << endl);

	// assign common vertices needed for construction.
	//TODO renumber if everything works
	const v v2(a, 0, 0);
	const v v3(a / 2, g, 0);
	const v v4(0, 0, h);
	const v v5(a, 0, h);
	const v v6(a / 2, g, h);
	const v v7(a / 2, g + s, 0);
	const v v8(-a / 2, g, 0);
	const v v9(-a / 2, g, h);
	const v v10(-a / 2, g + s, 0);
	const v v11(a / 2, g + s / 2, 0);
	const v v12(-a / 2, g + s / 2, 0);
	const v v13(a / 2 + b, g + s / 2, 0);

	const v v14(a / 2, g + s, h);
	const v v15(-a / 2, g + s, h);
	const v v16(-a / 2, g + s / 2, h);
	const v v17(a / 2, g + s / 2, h);

	const v v20(a / 2 + b, g + s / 2, h);
	const v v21(a / 2 + b, g + s / 2, 0);
	const v v22(a / 2 + b, g + s, h);

	///// intialisation of top/bottom segments
	obenInnen << origin << v2 << v3 << v4 << v5 << v6;
	obenAussenPrism << v3 << v6 << v7 << v8 << v9 << v10;
	obenAussenPr_rightTetrahedron << v3 << v11 << v13 << v6;
	obenAussenPr2T_prism << v3 << v6 << v11 << v8 << v9 << v12;

	// mirror right tetrahedron at x-axis to obtain left tetrahedron
	obenAussenPr_leftTetrahedrson = mirror(obenAussenPr_rightTetrahedron,
			xAxis);

	///// intialisiation of middle segments
	mitteAussenHexahedron << v3 << v6 << v14 << v11 << v8 << v9 << v15 << v12;
	mitteAussenP1 << v9 << v16 << v10 << v6 << v17 << v7;
	mitteAussenP2 << v8 << v9 << v10 << v3 << v6 << v7;

	/// this is tetrahedron of mitteAussenH2Pr right of mitteAussenHexahedron
	mitteAussenH2Pr_tetrahedron << v6 << v17 << v20 << v7;
	mitteAussen2PrH_tetrahedron = mirror(mitteAussenH2Pr_tetrahedron, xAxis);

	mitteAussenH2Pr_pyramid << v3 << v6 << v14 << v11 << v21;
	mitteAussen2PrH_pyramid = mirror(mitteAussenH2Pr_pyramid, xAxis);
}

number Generator::getVolume() const {
	return 3 * sqrt(3) / 2 * a * a * h + 3 * a * h * s
			+ 3 * sqrt(3) / 4 * s * s * h;
}

number Generator::getSurface() const {
	return 3 * a * a * sqrt(3) + 6 * a * sqrt(1 / 9 * h * h + s * s)
			+ 6 * sqrt(1 / 9 * h * h + s * s / 4) * (2 * a + s * sqrt(3));
}

void Generator::createGeometricObject(const CoordsArray & posIn) {
	uint numCoords = posIn.size();
	UG_ASSERT(
			numCoords == Tetrahedron || numCoords == Pyramid || numCoords == Prism || numCoords == Hexahedron,
			"wrong number of coordinates in posIn");

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
