/*
 * generator.cpp
 *
 *  Created on: 08.12.2011
 *      Author: marscher
 */
#include "generator.h"
namespace tkdGenerator {

TKDGeometryGenerator::TKDGeometryGenerator(number height, number baseEdgeLength,
		number diameter, number d_lipid) :
		h_corneocyte(height / 3), a_corneocyte(baseEdgeLength), w_corneocyte(
				diameter), d_lipid(d_lipid), R(0), index(0) {
	// TODO reserve memory
	indsOut.reserve(405);
	posOut.reserve(342);

	number a1 = sqrt(
			1.0 / 9.0 * height * height
					+ 1.0 / 3.0 * pow((diameter - 2.0 * baseEdgeLength), 2));

	number alpha = acos((diameter - 2.0 * baseEdgeLength) / (2.0 * a1));
	number beta = 90.0 / 180.0 * M_PI
			+ acos(1.0 / 3.0 * height / (a1 * sin(alpha)));
	number gamma = acos(1.0 / 3.0 * height / a1) + 90.0 / 180.0 * M_PI;

	number m1 = (d_lipid / 2) / tan(beta / 2);
	number m2 = (d_lipid / 2) / tan(gamma / 2);

	a_lipid = (sqrt(3) + a_corneocyte + m1 + m2) / sqrt(3);
	h_lipid = (h_corneocyte + d_lipid) / 3.0;

	UG_LOG(
			"h_c: " << h_corneocyte*3 << "\nd_l: " << d_lipid<<"\nh_l: " << h_lipid*3 << endl);

	//TODO is this correct?
	number w_lipid = w_corneocyte + d_lipid;

	s_lipid = 1 / sqrt(3) * (w_lipid - 2 * a_lipid);
	s_corneocyte = 1 / sqrt(3) * (w_corneocyte - 2 * a_corneocyte);

	initGeometricParams();
}

/**
 * creates top and bottom part of inner tkd with given offset and rotationOffset
 * @param offset
 * @param rotationOffset
 */
void TKDGeometryGenerator::createCorneocyteTop(const vector3& offset,
		const number rotationOffset) {
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

void TKDGeometryGenerator::createLipidTop(const vector3& offset,
		const number rotationOffset) {
	myTransform t(R, offset);

	// create G(Ki -> ObenInnen) = 6 prism with equilateral sites
	// step angle 60°
	for (int r = rotationOffset; r < 360 + rotationOffset; r += 60) {
		// sets the angle of rotation matrix which is applied on all points in this loop
		R.setAngle(r);
		createGeometricObject(t.perform(lipidTop));
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

void TKDGeometryGenerator::createLipidMiddle(const vector3& offset,
		const number rotationOffset) {
	myTransform t(R, offset);

	for (int r = rotationOffset; r < 360; r += 120) {
		R.setAngle(r);
		createGeometricObject(t.perform(sideQuadHexahedron));
	}
}

/**
 * creates two nested tetrakaidecahedrons (inner one with instanced geometric parameters)
 * and outer one with distance d_lipid / 2 to inner one. Starting construction in (0, 0, h/2),
 * so that center of geometry of tkd is (0, 0, 0)
 */
void TKDGeometryGenerator::createDomain() {
	vector3 offset(0, 0, h_corneocyte / 2);
	createCorneocyte(offset);
	createLipid(offset);
}

/**
 * creates a tkd starting construction with offset. Construction is performed down the z axis
 */
void TKDGeometryGenerator::createCorneocyte(const vector3& offset) {
	createCorneocyteTop(offset, 0);

	vector3 offset_h(offset.x, offset.y, offset.z - h_corneocyte);

	createCorneocyteMiddle(offset_h);

	// offset is same like middle because bottom is mirrored around x axis
	R.setMirrorZAxis(true);
	// in fact this is the bottom part, rotated with 60° relative to top
	createCorneocyteTop(offset_h, 60);
	// reset mirror status
	R.setMirrorZAxis(false);
}

void TKDGeometryGenerator::createLipid(const vector3& offset) {
	createLipidTop(vector3(offset.x, offset.y, offset.z + h_corneocyte));

	// create middle
	// todo fix offset
	createLipidMiddle(vector3(offset.x, offset.y, offset.z + h_corneocyte));
	R.setMirrorZAxis(true);
	// same elements but rotated by 60° and mirrored at z axis
	createLipidMiddle(
			vector3(offset.x, offset.y, offset.z - 2 * h_corneocyte / 2), 60);

	// bottom
	createLipidTop(vector3(offset.x, offset.y, offset.z - 2 * h_corneocyte));
	// reset mirror status
	R.setMirrorZAxis(false);
}

/**
 * init basis geometrical parameters
 */
void TKDGeometryGenerator::initGeometricParams() {
	// height of base triangle of top inner prism
	number g_cornoecyte = sqrt(3) * a_corneocyte / 2;
	number g_lipid = sqrt(3) * a_lipid / 2;
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

	//// lipid vertices
	// top
	const vector3 v19(0, 0, d_lipid / 2);
	const vector3 v20(a_lipid, 0, v19.z);
	const vector3 v21(a_lipid / 2, g_lipid, v19.z);
	// side hexagon hexahedron
	//TODO find correct values
	number x = 0;
	const vector3 v22(a_lipid / 2, g_lipid + x, h_lipid);
	const vector3 v23(a_lipid, v22.y, 0);
	const vector3 v24(a_lipid, 0, v19.z);

	///// Initialization of top/bottom segments
	obenInnen << origin << v1 << v4 << v2 << v3 << v5;
	obenAussenPrism << v4 << v5 << v8 << v6 << v7 << v9;
	obenAussenPr_rightTetrahedron << v4 << v12 << v16 << v5;
	obenAussenPr2T_prism << v4 << v5 << v12 << v6 << v7 << v13;

	// mirror right tetrahedron at x-axis to obtain left tetrahedron
	obenAussenPr_leftTetrahedrson = mirror(obenAussenPr_rightTetrahedron,
			xAxis);

	///// Initialization of middle segments
	mitteAussenHexahedron << v4 << v5 << v10 << v12 << v6 << v7 << v11 << v13;
	mitteAussenP1 << v7 << v14 << v9 << v5 << v15 << v8;
	mitteAussenP2 << v6 << v7 << v9 << v4 << v5 << v8;

	/// this is tetrahedron of mitteAussenH2Pr right of mitteAussenHexahedron
	mitteAussenH2Pr_tetrahedron << v5 << v15 << v17 << v8;
	mitteAussen2PrH_tetrahedron = mirror(mitteAussenH2Pr_tetrahedron, xAxis);

	mitteAussenH2Pr_pyramid << v4 << v5 << v10 << v12 << v16;
	mitteAussen2PrH_pyramid = mirror(mitteAussenH2Pr_pyramid, xAxis);

	///// Initialization of lipid
	/// top segment
	lipidTop << origin << v1 << v4 << v19 << v20 << v21;
	//	sideHexagonHexahedron;

	// shares face with ObenAussenPrism
	sideQuadHexahedron << v5 << v8 << v9 << v7;
}

number TKDGeometryGenerator::getVolume(int subset) const {
	if (subset == CORNEOCYTE)
		return getVolume(a_corneocyte, s_corneocyte, h_corneocyte);
	if (subset == LIPID)
		return getVolume(a_lipid, s_lipid, h_lipid);
	return -1;
}

number TKDGeometryGenerator::getSurface(int subset) const {
	if (subset == CORNEOCYTE)
		return getSurface(a_corneocyte, s_corneocyte, h_corneocyte);
	if (subset == LIPID)
		return getSurface(a_lipid, s_lipid, h_lipid);
	return -1;
}

number TKDGeometryGenerator::getVolume(number a, number s, number h) const {
	return 3 * sqrt(3) / 2 * a * a * h + 3 * a * h * s
			+ 3 * sqrt(3) / 4 * s * s * h;
}

number TKDGeometryGenerator::getSurface(number a, number s, number h) const {
	return 3 * a * a * sqrt(3) + 6 * a * sqrt(1.0 / 9 * h * h + s * s)
			+ 6 * sqrt(1.0 / 9 * h * h + s * s / 4) * (2 * a + s * sqrt(3));
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
