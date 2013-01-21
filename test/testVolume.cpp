/*
 * testVolume.cpp
 *
 *  Created on: 12.01.2012
 *      Author: marscher
 */
#include <boost/test/unit_test.hpp>
#include <vector>
#include "lib_grid/lib_grid.h"
#include "../util/volume_calculation.h"

using namespace ug;
//using std::vector;
using namespace std;

/**
 * common code for every volume testsuite testcase.
 * It sets up a Grid instance and a VertexAttachmentAccessor
 */
struct volFixture {
	volFixture() :
			grid(VRTOPT_STORE_ASSOCIATED_FACES) {
		grid.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(grid, aPosition);
	}
	Grid grid;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	VolumeDescriptor vd;
	vector<VertexBase*> vertices;
};

BOOST_FIXTURE_TEST_SUITE(volume, volFixture);

BOOST_AUTO_TEST_CASE(volHexahedron) {
	for(uint i = 0; i < 8; i++)
	vertices.push_back(*grid.create<Vertex>());

	//	assign coordinates of unit cube
	aaPos[vertices[0]] = vector3(1, 0, 0);
	aaPos[vertices[1]] = vector3(1, 1, 0);
	aaPos[vertices[2]] = vector3(-1, 1, 0);
	aaPos[vertices[3]] = vector3(-1, 0, 0);

	aaPos[vertices[4]] = vector3(1, 0, 1);
	aaPos[vertices[5]] = vector3(1, 1, 1);
	aaPos[vertices[6]] = vector3(-1, 1, 1);
	aaPos[vertices[7]] = vector3(-1, 0, 1);

	vd.set_num_vertices(vertices.size());
	for (uint i = 0; i < vertices.size(); ++i)
	vd.set_vertex(i, vertices[i]);

	//	create the hexahedron
	Hexahedron* h = *grid.create<Hexahedron>(vd);
	// calculate its volume
	number vol = CalculateVolume(*h, aaPos);
	// check if calculated volume equals 1
	BOOST_CHECK_CLOSE(vol, 2.0, 0.001);
}

BOOST_AUTO_TEST_CASE(volTetrahedron) {
	//	create vertices
	for(uint i = 0; i < 4; i++)
	vertices.push_back(*grid.create<Vertex>());

	//	assign coordinates
	aaPos[vertices[0]] = vector3(0, 1, 0);
	aaPos[vertices[1]] = vector3(1, 0, 0);
	aaPos[vertices[2]] = vector3(-1, 0, 0);
	aaPos[vertices[3]] = vector3(0, 0.5, 1);

	vd.set_num_vertices(vertices.size());
	for (uint i = 0; i < vertices.size(); ++i)
	vd.set_vertex(i, vertices[i]);

	//	create the tetrahedron
	Tetrahedron* t = *grid.create<Tetrahedron>(vd);
	number vol = CalculateVolume(*t, aaPos);
	UG_LOG("vol of tetrahedron: " << vol << endl)
	BOOST_CHECK_CLOSE(vol, 1./3, 0.001);
}

BOOST_AUTO_TEST_CASE(volPrism) {
	//	create vertices
	for(uint i = 0; i < 6; i++)
		vertices.push_back(*grid.create<Vertex>());

	//	assign coordinates
	aaPos[vertices[0]] = vector3(0, 0, 0);
	aaPos[vertices[1]] = vector3(1, 0, 0);
	aaPos[vertices[2]] = vector3(0, 1, 0);

	aaPos[vertices[3]] = vector3(0, 0, 1);
	aaPos[vertices[4]] = vector3(1, 0, 1);
	aaPos[vertices[5]] = vector3(0, 1, 1);

	vd.set_num_vertices(vertices.size());
	for (uint i = 0; i < vertices.size(); ++i)
		vd.set_vertex(i, vertices[i]);

	//	create the prism
	Prism* p = *grid.create<Prism>(vd);

	number vol = CalculateVolume(*p, aaPos);
	UG_LOG("vol of prism: " << vol << endl);

	BOOST_CHECK_CLOSE(vol, 0.5, 0.001);
}

BOOST_AUTO_TEST_CASE(volPyramid) {
	//	create vertices
	for(uint i = 0; i < 5; i++)
	vertices.push_back(*grid.create<Vertex>());

	//	assign coordinates
	aaPos[vertices[0]] = vector3(0, 0, 0);
	aaPos[vertices[1]] = vector3(1, 0, 0);
	aaPos[vertices[2]] = vector3(1, 1, 0);
	aaPos[vertices[3]] = vector3(0, 1, 0);
	// top
	aaPos[vertices[4]] = vector3(0.5, 0.5, 1);

	vd.set_num_vertices(vertices.size());
	for (uint i = 0; i < vertices.size(); ++i)
	vd.set_vertex(i, vertices[i]);

	//	create the prism
	Pyramid* p = *grid.create<Pyramid>(vd);

	number vol = CalculateVolume(*p, aaPos);
	UG_LOG("vol of pyramid: " << vol << endl);
	BOOST_CHECK_CLOSE(vol, 1./3, 0.001);

	for(uint i = 0; i < 5; i++) {
		VecScale(aaPos[vertices[i]],aaPos[vertices[i]], 2.0);
	}
	vol = CalculateVolume(*p, aaPos);
	BOOST_CHECK_CLOSE(vol, 8.0 / 3.0, 0.001);
}

BOOST_AUTO_TEST_SUITE_END();
