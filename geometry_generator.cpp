/*
 * generator.cpp
 *
 *  Created on: 08.12.2011
 *      Author: marscher
 */
#include "geometry_generator.h"
#include "util/geometric_helper.h"
#include "coordinates.h"
#include <algorithm>

namespace ug {
namespace tkd {

const vector3 TKDGeometryGenerator::origin(0, 0, 0);

TKDGeometryGenerator::TKDGeometryGenerator() :
		b_createLipid(true), R(0), m_num_verts_created(0) {
}

TKDGeometryGenerator::TKDGeometryGenerator(number a, number w, number h,
		bool createLipid, number d_lipid) :
		b_createLipid(createLipid), R(0), m_num_verts_created(0) {

	setGeometricParams(a, w, h, d_lipid);
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
 * creates two nested tkd starting construction in (0, 0, h/2),
 * so that center of geometry is (0, 0, 0)
 */
void TKDGeometryGenerator::createGeometry() {
	// clear output data
	m_mCoordIndexMapping.clear();
	indsOut.clear();

	// init basic volume coordinates
	initGeometricParams();

	vector3 offset(0, 0, h_corneocyte / 2);
	createCorneocyte(offset);

	if(b_createLipid) {
		// create top part. Due to stolen coordinate calculation,
		// offset starts with 60° instead of 0.
		createLipidMatrix(origin, 60);
		// and bottom part
		createLipidMatrix(origin, 120, true);

		// clear lipid arrays after creation
		obenInnenPrismL.clear();
		obenAussen_leftPrismL.clear();
		obenAussen_rightPrismL.clear();
		upperHexahedronL.clear();
		bottomLeftPrismL.clear();
		bottomOuterHexahedronL.clear();
		bottomRightPrismL.clear();
		sideQuad.clear();
	}
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

	// clear memory of corn arrays
	obenInnen.clear();
	obenAussenPrism.clear();
	obenAussenPr_rightTetrahedron.clear();
	obenAussenPr_leftTetrahedron.clear();
	obenAussenPr2T_prism.clear();
	mitteAussenP1.clear();
	mitteAussenP2.clear();
	mitteAussenH2Pr_tetrahedron.clear();
	mitteAussenH2Pr_pyramid.clear();
	mitteAussen2PrH_tetrahedron.clear();
	mitteAussen2PrH_pyramid.clear();
	mitteAussenHexahedron.clear();
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
	m_num_verts_created = 0;
	indsOut.reserve(114);

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
	obenInnen << CLEAR << origin << v1 << v4 << v2 << v3 << v5;
	obenAussenPrism << CLEAR << v4 << v5 << v8 << v6 << v7 << v9;
	obenAussenPr_rightTetrahedron << CLEAR << v16 << v12 << v4 << v5;
	obenAussenPr2T_prism << CLEAR << v4 << v5 << v12 << v6 << v7 << v13;

	// mirror right tetrahedron at x-axis to obtain left tetrahedron
	obenAussenPr_leftTetrahedron = mirror(obenAussenPr_rightTetrahedron, xAxis);
	// flip orientation after mirroring
	flipOrientationTetrahedron(obenAussenPr_leftTetrahedron);

	///// Initialization of middle segments
	mitteAussenHexahedron << CLEAR << v5 << v7 << v11 << v10 << v4 << v6 << v13
			<< v12;

	mitteAussenP1 << CLEAR << v5 << v15 << v8 << v7 << v14 << v9;
	mitteAussenP2 << CLEAR << v4 << v5 << v8 << v6 << v7 << v9;

	/// this is tetrahedron of mitteAussenH2Pr right of mitteAussenHexahedron
	mitteAussenH2Pr_tetrahedron << CLEAR << v5 << v15 << v17 << v8;
	mitteAussen2PrH_tetrahedron = mirror(mitteAussenH2Pr_tetrahedron, xAxis);
	// flip orientation after mirroring
	flipOrientationTetrahedron(mitteAussen2PrH_tetrahedron);
	mitteAussenH2Pr_pyramid << CLEAR << v10 << v5 << v4 << v12 << v16;
	mitteAussen2PrH_pyramid = mirror(mitteAussenH2Pr_pyramid, xAxis);
	// fix orientation
	flipOrientationPyramid(mitteAussen2PrH_pyramid);

	///// lipid matrix
	if(b_createLipid)
		initLipidGeometric();
}

void TKDGeometryGenerator::initLipidGeometric() {
	CoordsArray c(76);
	CalculateLipidCoords(c, a_corneocyte, h_corneocyte * 3, w_corneocyte,
			d_lipid, vector3(0, 0, 3 * h_corneocyte / 2));

	// inner point 7, 8, 13
	obenInnenPrismL << CLEAR << c[13] << c[7] << c[8] << c[6] << c[0] << c[1];
	// inner point 26, 27, 7, 8
	sideQuad << CLEAR << c[26] << c[27] << c[8] << c[7] << c[14] << c[15]
			<< c[1] << c[0];

	// inner points 32, 10, 31
	obenAussen_leftPrismL << CLEAR << c[32] << c[10] << c[31] << c[20] << c[3]
			<< c[19];
	obenAussen_rightPrismL = mirror(obenAussen_leftPrismL, xAxis);
	flipOrientationPrism(obenAussen_rightPrismL);

	// inner points 32, 31, 45
	bottomLeftPrismL << CLEAR << c[32] << c[31] << c[45] << c[20] << c[19]
			<< c[57];
	bottomRightPrismL = mirror(bottomLeftPrismL, xAxis);
	flipOrientationPrism(bottomRightPrismL);

	// inner points 27, 26, 39, 40
	bottomOuterHexahedronL << CLEAR << c[27] << c[26] << c[39] << c[40] << c[15]
			<< c[14] << c[51] << c[52];
	// inner points 10, 11, 32, 33
	upperHexahedronL << CLEAR << c[11] << c[10] << c[32] << c[33] << c[4]
			<< c[3] << c[20] << c[21];
}

number TKDGeometryGenerator::getVolume(int subset) const {
	if (subset == CORNEOCYTE)
		return getVolume(a_corneocyte, s_corneocyte, 3.0 * h_corneocyte);
	else if (subset == LIPID) {
		if(b_createLipid)
			return getVolume(a_lipid, s_lipid, 3.0 * h_lipid)
					- getVolume(CORNEOCYTE);
		else
			return -1;
	} else
		UG_THROW("no valid subset given.");
}

number TKDGeometryGenerator::getSurface(int subset) const {
	if (subset == CORNEOCYTE || subset == BOUNDARY_CORN)
		return getSurface(a_corneocyte, s_corneocyte, 3 * h_corneocyte);
	else if (subset == LIPID || subset == BOUNDARY_LIPID)
		if(b_createLipid)
			return getSurface(a_lipid, s_lipid, 3 * h_lipid);
		else
			return -1;
	else
		UG_THROW("no valid subset given.");
}

number TKDGeometryGenerator::getVolume(number a, number s, number h) {
	number sq = sqrt(3) * a + s;
	return 0.5 * sqrt(3) * h * sq * sq;
}

number TKDGeometryGenerator::getSurface(number a, number s, number h) {
	return 3.0 * sqrt(3.0) * a * a + 6.0 * a * sqrt(1.0 / 9.0 * h * h + s * s)
			+ 6.0 * (2.0 * a + sqrt(3.0) * s)
					* sqrt(1.0 / 9.0 * h * h + 1.0 / 4.0 * s * s);
}

TKDGeometryGenerator::CoordIndexMap&
TKDGeometryGenerator::getPositions() {
	return m_mCoordIndexMapping;
}

number TKDGeometryGenerator::getHeight() const {
	return h_lipid * 3.0;
}

number TKDGeometryGenerator::getOverlap() const {
	return s_lipid;
}

void TKDGeometryGenerator::createGeometricObject(const CoordsArray& posIn) {
	uint numCoords = posIn.size();
	UG_ASSERT(numCoords == Tetrahedron || numCoords == Pyramid ||
			numCoords == Prism || numCoords == Hexahedron,
			"wrong number of coordinates in posIn");

	// how much indices are belonging to this geometric object
	indsOut.push_back(numCoords);

	// lookup index by position i, if position not yet known.
	int index = -1;
	for (uint i = 0; i < numCoords; ++i) {
		const vector3& pos = posIn[i];

		// try to find the position
		CoordIndexMap::right_const_iterator li =
				m_mCoordIndexMapping.right.find(posIn[i]);
		// position not yet known, create new index and insert into map
		if (li == m_mCoordIndexMapping.right.end()) {
			index = m_num_verts_created++;
			m_mCoordIndexMapping.insert(CoordIndexMap::value_type(index, pos));
		} else {
			// lookup index in bimap per position;
			index = m_mCoordIndexMapping.right.at(pos);
		}

		indsOut.push_back(index);
	}
}

void TKDGeometryGenerator::flipOrientationPrism(CoordsArray& prismPos) const {
	UG_ASSERT(prismPos.size() == Prism, "no prism given.");
	std::swap(prismPos[0], prismPos[3]);
	std::swap(prismPos[1], prismPos[4]);
	std::swap(prismPos[2], prismPos[5]);
}

void TKDGeometryGenerator::flipOrientationTetrahedron(
		CoordsArray& tetrahedronPos) const {
	UG_ASSERT(tetrahedronPos.size() == Tetrahedron, "no tetrahedron given.");
	std::swap(tetrahedronPos[0], tetrahedronPos[2]);
}

void TKDGeometryGenerator::flipOrientationPyramid(CoordsArray& pyramidPos) const {
	UG_ASSERT(pyramidPos.size() == Pyramid, "no pyramid given.");
	std::swap(pyramidPos[0], pyramidPos[3]);
	std::swap(pyramidPos[1], pyramidPos[2]);
}

void TKDGeometryGenerator::flipOrientationHexahedron(CoordsArray& hexaPos) const {
	UG_ASSERT(hexaPos.size() == Hexahedron, "no hexahedron given.");
	std::swap(hexaPos[0], hexaPos[4]);
	std::swap(hexaPos[1], hexaPos[5]);
	std::swap(hexaPos[2], hexaPos[6]);
	std::swap(hexaPos[3], hexaPos[7]);
}

void TKDGeometryGenerator::setHeight(number height) {
	// we use one third of h internally
	this->h_corneocyte = height / 3.0;
}

void TKDGeometryGenerator::setBaseEdgeLength(number a) {
	this->a_corneocyte = a;
}

void TKDGeometryGenerator::setDiameter(number w) {
	this->w_corneocyte = w;
}

void TKDGeometryGenerator::setLipidThickness(number d) {
	this->d_lipid = d;
}

void TKDGeometryGenerator::setGeometricParams(number a, number w, number h,
		number d_lipid) {
	setBaseEdgeLength(a);
	setHeight(h);
	setDiameter(w);
	setLipidThickness(d_lipid);

	// update overlap as it depends on a and w
	updateOverlap(CORNEOCYTE);
	if(b_createLipid)
		setLipidParameters();
}

void TKDGeometryGenerator::updateOverlap(int subset) {
	if (subset == CORNEOCYTE) {
		this->s_corneocyte = (sqrt(3.0) / 3.0)
				* (w_corneocyte - 2 * a_corneocyte);
		UG_DLOG(APP, 1,
				"updated corneocyte overlap to: " << s_corneocyte << std::endl);
	} else if (subset == LIPID) {
		this->s_lipid = (sqrt(3.0) / 3.0) * (w_lipid - 2.0 * a_lipid);
		UG_DLOG(APP, 1,
				"updated lipid overlap to: " << s_lipid << std::endl);
	}
}

void TKDGeometryGenerator::setLipidBaseEdgeLength() {
	number h = 3 * h_corneocyte;
	// 1/9 h^2 + s^2
	number a1 = sqrt(1 / 9. * h * h + s_corneocyte * s_corneocyte);

	number alpha = acos((w_corneocyte - 2 * a_corneocyte) / (2 * a1));
	number beta = 1.0 / 2.0 * M_PI + acos(h / (3.0 * a1 * sin(alpha)));
	number gamma = 1.0 / 2.0 * M_PI + acos(h / (3.0 * a1));

	UG_DLOG(APP, 1,
			"alpha: " << alpha << " beta: " << beta << " gamma: " << gamma << std::endl)

	number m1 = d_lipid / (2.0 * tan(beta / 2.0));
	number m2 = d_lipid / (2.0 * tan(gamma / 2.0));

	this->a_lipid = sqrt(3) / 3.0 * (sqrt(3) * a_corneocyte + m1 + m2);
}

void TKDGeometryGenerator::setLipidHeight() {
	this->h_lipid = (3.0 * h_corneocyte + d_lipid) / 3.0;
}

void TKDGeometryGenerator::setLipidDiameter() {
	this->w_lipid = (h_lipid * 3.0) * (w_corneocyte - 2.0 * a_corneocyte)
			/ (h_corneocyte * 3.0) + 2.0 * a_lipid;
}

void TKDGeometryGenerator::setLipidParameters() {
	setLipidBaseEdgeLength();
	setLipidHeight();
	setLipidDiameter();

	updateOverlap (LIPID);

	UG_DLOG(APP, 1,
			"a_l: " << a_lipid << " w_l: " << w_lipid << " s_l: " << s_lipid << std::endl);
}

bool TKDGeometryGenerator::createLipid() const {
	return b_createLipid;
}

void TKDGeometryGenerator::setCreateLipid(bool createLipid) {
	this->b_createLipid = createLipid;
	if(!this->b_createLipid) {
		initLipidGeometric();
	}
}

}// end of namespace tkd
}// end of namespace ug
