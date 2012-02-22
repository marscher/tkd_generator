/*
 * generator.cpp
 *
 *  Created on: 08.12.2011
 *      Author: marscher
 */
#include "generator.h"
#include "geometric_helper.h"
#include "coordinates.h"
#include <stdexcept>
#include <algorithm>

namespace tkdGenerator {
using namespace std;

TKDGeometryGenerator::TKDGeometryGenerator(number height, number baseEdgeLength,
		number diameter, number d_lipid) :
		h_corneocyte(height / 3), a_corneocyte(baseEdgeLength),
				w_corneocyte(diameter), d_lipid(d_lipid), R(0), index(0) {
	indsOut.reserve(702);
	posOut.reserve(819);

	number w2a = w_corneocyte - 2 * a_corneocyte;
	s_corneocyte = 1 / sqrt(3) * w2a;

	number a1 = sqrt(1 / 9. * height * height + 1 / 3. * w2a * w2a);

	number alpha = acos(w2a / (2 * a1));
	number beta = 90 / 180. * M_PI + acos(1 / 3. * height / (a1 * sin(alpha)));
	number gamma = acos(1 / 3. * height / a1) + 90 / 180. * M_PI;

	number m1 = (d_lipid / 2) / tan(beta / 2);
	number m2 = (d_lipid / 2) / tan(gamma / 2);

	a_lipid = (sqrt(3) + a_corneocyte + m1 + m2) / sqrt(3);
	// height of lipid is full height of inner tkd + d_lipid
	h_lipid = (height + d_lipid) / 3;

	//TODO is this correct?
	number w_lipid = w_corneocyte + d_lipid;

	s_lipid = 1 / sqrt(3) * (w_lipid - 2 * a_lipid);

	initGeometricParams();
	count = 0;
}

/**
 * creates top and bottom part of inner tkd with given offset and rotationOffset
 * @param offset
 * @param rotationOffset
 */
void TKDGeometryGenerator::createCorneocyteTop(const vector3& offset,
		const number rotationOffset, bool bottom) {

	// if we are creating the bottom part, flip orientation of all volumes
	if (bottom) {
		flipOrientationPrism(obenInnen);
		flipOrientationPrism(obenAussenPrism);
		flipOrientationPrism(obenAussenPr2T_prism);

		flipOrientationTetrahedron(obenAussenPr_leftTetrahedron);
		flipOrientationTetrahedron(obenAussenPr_rightTetrahedron);
		// mirror at z axis
		R.setMirrorZAxis(true);
	}
	// transform all geometric objects with R and offset.
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
		createGeometricObject(t.perform(obenAussenPr_leftTetrahedron));

		/* tetrahedron right of prism ObenAussenPr2T which shares
		 3 points with prism obenAussenPr2T_prism */
		createGeometricObject(t.perform(obenAussenPr_rightTetrahedron));
	}

	// reset mirror status
	if (bottom) {
		R.setMirrorZAxis(false);
	}
}

void TKDGeometryGenerator::createCorneocyteMiddle(const vector3& origin) {
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
 * creates one tkd starting construction in (0, 0, h/2),
 * so that center of geometry of tkd is (0, 0, 0)
 */
void TKDGeometryGenerator::createDomain() {
	vector3 offset(0, 0, h_corneocyte / 2);
	createCorneocyte(offset);

	RotationMatrix R(60);
	for (uint i = 0; i < posOut.size(); i++) {
		posOut[i] = R * posOut[i];
	}

	createLipidMatrix(origin);
	createLipidMatrix(origin, 60, true);
}

/**
 * creates a tkd starting construction with offset. Construction is performed down the z axis
 */
void TKDGeometryGenerator::createCorneocyte(const vector3& offset) {
	createCorneocyteTop(offset, 0);

	vector3 offset_h(offset.x, offset.y, offset.z - h_corneocyte);

	createCorneocyteMiddle(offset_h);

	// offset is same like middle because bottom is mirrored around z axis
	// in fact this is the bottom part, rotated with 60° relative to top
	createCorneocyteTop(offset_h, 60, true);
}

void TKDGeometryGenerator::createLipidMatrix(const vector3& offset,
		const number rotationOffset, bool bottom) {

	if (bottom) {
		flipOrientationPrism(obenInnenPrismL);
		flipOrientationPrism(obenAussen_leftPrismL);
		flipOrientationPrism(obenAussen_rightPrismL);
		flipOrientationHexahedron(upperHexahedronL);
		flipOrientationPrism(bottomLeftPrismL);
		flipOrientationHexahedron(bottomOuterHexahedronL);
		flipOrientationPrism(bottomRightPrismL);
		flipOrientationHexahedron(sideQuad);

		R.setMirrorZAxis(true);
	}

	myTransform t(R, offset);

	for (int r = rotationOffset; r < 360 + rotationOffset; r += 60) {
		R.setAngle(r);
		createGeometricObject(t.perform(obenInnenPrismL));
	}

	for (int r = rotationOffset; r < 360 + rotationOffset; r += 120) {
		R.setAngle(r);
		createGeometricObject(t.perform(sideQuad));
		createGeometricObject(t.perform(obenAussen_leftPrismL));
		createGeometricObject(t.perform(upperHexahedronL));
		createGeometricObject(t.perform(obenAussen_rightPrismL));
		createGeometricObject(t.perform(bottomLeftPrismL));
		createGeometricObject(t.perform(bottomOuterHexahedronL));
		createGeometricObject(t.perform(bottomRightPrismL));
	}

	if (bottom) {
		R.setMirrorZAxis(false);
	}
}

/**
 * init basis geometrical parameters
 */
void TKDGeometryGenerator::initGeometricParams() {
	// height of base triangle of top inner prism
	number g_cornoecyte = sqrt(3) * a_corneocyte / 2;
	// height of base triangle of tetrahedron of ObenAussenPr2T
	number b = sqrt(3) * s_corneocyte / 2;

	// assign common vertices needed for construction.
	const vector3 v1(a_corneocyte, 0, 0);
	const vector3 v2(0, 0, h_corneocyte);
	const vector3 v3(v1.x, 0, h_corneocyte);

	const vector3 v4(v1.x / 2, g_cornoecyte, 0);
	const vector3 v5(v4.x, v4.y, h_corneocyte);

	const vector3 v6(-v4.x, v4.y, 0);
	const vector3 v7(v6.x, v4.y, h_corneocyte);

	const vector3 v8(v4.x, v4.y + s_corneocyte, 0);
	const vector3 v9(v6.x, v8.y, 0);
	const vector3 v10(v4.x, v8.y, h_corneocyte);
	const vector3 v11(v6.x, v8.y, h_corneocyte);

	const vector3 v12(v4.x, v4.y + s_corneocyte / 2, 0);
	const vector3 v13(v6.x, v12.y, 0);
	const vector3 v14(v6.x, v12.y, h_corneocyte);
	const vector3 v15(v4.x, v12.y, h_corneocyte);

	const vector3 v16(v4.x + b, v12.y, 0);
	const vector3 v17(v16.x, v12.y, h_corneocyte);
	const vector3 v18(v16.x, v8.y, h_corneocyte);

	///// lipid vertices
	// top/bottom prism
	const vector3 v20(a_lipid, 0, d_lipid / 2 + h_corneocyte);
	const vector3 v21(a_lipid / 2.0, a_lipid * sqrt(3) / 2, v20.z);
	const vector3 v22(0, 0, v20.z);

	///// Initialization of top/bottom segments
	obenInnen << origin << v1 << v4 << v2 << v3 << v5;
	obenAussenPrism << v4 << v5 << v8 << v6 << v7 << v9;
	obenAussenPr_rightTetrahedron << v16 << v12 << v4 << v5;
	obenAussenPr2T_prism << v4 << v5 << v12 << v6 << v7 << v13;

	// mirror right tetrahedron at x-axis to obtain left tetrahedron
	obenAussenPr_leftTetrahedron = mirror(obenAussenPr_rightTetrahedron, xAxis);
	// flip orientation after mirroring
	flipOrientationTetrahedron(obenAussenPr_leftTetrahedron);

	///// Initialization of middle segments
	mitteAussenHexahedron << v5 << v7 << v11 << v10 << v4 << v6 << v13 << v12;

	mitteAussenP1 << v5 << v15 << v8 << v7 << v14 << v9;
	mitteAussenP2 << v4 << v5 << v8 << v6 << v7 << v9;

	/// this is tetrahedron of mitteAussenH2Pr right of mitteAussenHexahedron
	mitteAussenH2Pr_tetrahedron << v5 << v15 << v17 << v8;
	mitteAussen2PrH_tetrahedron = mirror(mitteAussenH2Pr_tetrahedron, xAxis);
	// flip orientation after mirroring
	flipOrientationTetrahedron(mitteAussen2PrH_tetrahedron);
	mitteAussenH2Pr_pyramid << v10 << v5 << v4 << v12 << v16;
	mitteAussen2PrH_pyramid = mirror(mitteAussenH2Pr_pyramid, xAxis);
	// fix orientation
	flipOrientationPyramid(mitteAussen2PrH_pyramid);

	///// lipid matrix
	// fixme substitute inner points with v_i
	CoordsArray c = CalculateLipidCoords(a_corneocyte, h_corneocyte * 3,
			w_corneocyte, d_lipid, vector3(0, 0, 3 * h_corneocyte / 2));

	// inner point 7, 8, 13
	obenInnenPrismL << c[13] << c[7] << c[8] << c[6] << c[0] << c[1];
	// inner point 26, 27, 7, 8
	sideQuad << c[26] << c[27] << c[8] << c[7] << c[14] << c[15] << c[1]
			<< c[0];

	// inner points 32, 10, 31
	obenAussen_leftPrismL << c[32] << c[10] << c[31] << c[20] << c[3] << c[19];
	obenAussen_rightPrismL = mirror(obenAussen_leftPrismL, xAxis);
	flipOrientationPrism(obenAussen_rightPrismL);

	// inner points 32, 31, 45
	bottomLeftPrismL << c[32] << c[31] << c[45] << c[20] << c[19] << c[57];
	bottomRightPrismL = mirror(bottomLeftPrismL, xAxis);
	flipOrientationPrism(bottomRightPrismL);

	// inner points 27, 26, 39, 40
	bottomOuterHexahedronL << c[27] << c[26] << c[39] << c[40] << c[15] << c[14]
			<< c[51] << c[52];
	// inner points 10, 11, 32, 33
	upperHexahedronL << c[11] << c[10] << c[32] << c[33] << c[4] << c[3]
			<< c[20] << c[21];
}

number TKDGeometryGenerator::getVolume() const {
	return getVolume(a_corneocyte, s_corneocyte, h_corneocyte);
}

number TKDGeometryGenerator::getSurface() const {
	return getSurface(a_corneocyte, s_corneocyte, h_corneocyte);
}

number TKDGeometryGenerator::getVolume(number a, number s, number h) const {
	return 3 * sqrt(3) / 2 * a * a * h + 3 * a * h * s
			+ 3 * sqrt(3) / 4 * s * s * h;
}

number TKDGeometryGenerator::getSurface(number a, number s, number h) const {
	return 3 * a * a * sqrt(3) + 6 * a * sqrt(1. / 9 * h * h + s * s)
			+ 6 * sqrt(1. / 9 * h * h + s * s / 4) * (2 * a + s * sqrt(3));
}

const IndexArray& TKDGeometryGenerator::getIndices() const {
	return indsOut;
}

const CoordsArray& TKDGeometryGenerator::getPositions() const {
	return posOut;
}

number TKDGeometryGenerator::getHeight() const {
	return h_lipid * 3;
}

number TKDGeometryGenerator::getOverlap() const {
	return s_lipid;
}

void TKDGeometryGenerator::createGeometricObject(const CoordsArray& posIn) {
	count++;
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

void TKDGeometryGenerator::flipOrientationPrism(CoordsArray& prismPos) {
	UG_ASSERT(prismPos.size() == Prism, "no prism given.");
	swap(prismPos[0], prismPos[3]);
	swap(prismPos[1], prismPos[4]);
	swap(prismPos[2], prismPos[5]);
}

void TKDGeometryGenerator::flipOrientationTetrahedron(
		CoordsArray& tetrahedronPos) {
	UG_ASSERT(tetrahedronPos.size() == Tetrahedron, "no tetrahedron given.");
	swap(tetrahedronPos[0], tetrahedronPos[2]);
}

void TKDGeometryGenerator::flipOrientationPyramid(CoordsArray& pyramidPos) {
	UG_ASSERT(pyramidPos.size() == Pyramid, "no pyramid given.");
	swap(pyramidPos[0], pyramidPos[3]);
	swap(pyramidPos[1], pyramidPos[2]);
}

void TKDGeometryGenerator::flipOrientationHexahedron(CoordsArray& hexaPos) {
	UG_ASSERT(hexaPos.size() == Hexahedron, "no hexahedron given.");
	swap(hexaPos[0], hexaPos[4]);
	swap(hexaPos[1], hexaPos[5]);
	swap(hexaPos[2], hexaPos[6]);
	swap(hexaPos[3], hexaPos[7]);
}

} // end of namespace tkdGenerator
