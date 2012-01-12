/*
 * testVolume.cpp
 *
 *  Created on: 12.01.2012
 *      Author: marscher
 */

#include <boost/test/unit_test.hpp>
#include "lib_grid/lib_grid.h"
#include "../util/volume_calculation.h"
using namespace ug;

BOOST_AUTO_TEST_SUITE(volume);

BOOST_AUTO_TEST_CASE(volHexahedron) {
	Grid grid;

	//	create vertices
	Vertex* v1 = *grid.create<Vertex> ();
	Vertex* v2 = *grid.create<Vertex> ();
	Vertex* v3 = *grid.create<Vertex> ();
	Vertex* v4 = *grid.create<Vertex> ();
	Vertex* v5 = *grid.create<Vertex> ();
	Vertex* v6 = *grid.create<Vertex> ();
	Vertex* v7 = *grid.create<Vertex> ();
	Vertex* v8 = *grid.create<Vertex> ();

	//	attach position data
	grid.attach_to_vertices(aPosition);
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	//	assign coordinates of unit cube
	aaPos[v1] = vector3(0, -1, 0);
	aaPos[v2] = vector3(0, 1, 0);
	aaPos[v3] = vector3(1, 0, 0);
	aaPos[v4] = vector3(-1, 0, 0);

	aaPos[v5] = vector3(0, -1, 1);
	aaPos[v6] = vector3(0, 1, 1);
	aaPos[v7] = vector3(1, 0, 1);
	aaPos[v8] = vector3(-1, 0, 1);

	//	create the hexahedron
	HexahedronDescriptor hex(v1, v2, v3, v4, v5, v6, v7, v8);
	grid.create<Hexahedron> (hex);

	number vol = CalculateVolume(hex);
	BOOST_CHECK_EQUAL(vol, 1);
}

BOOST_AUTO_TEST_CASE(volTetrahedron) {
	BOOST_MESSAGE("testing vol calc of tetrahedron");
	Grid grid;
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid,aPosition);
	//	create vertices
	Vertex* v1 = *grid.create<Vertex> ();
	Vertex* v2 = *grid.create<Vertex> ();
	Vertex* v3 = *grid.create<Vertex> ();
	Vertex* v4 = *grid.create<Vertex> ();

	//	assign coordinates
	aaPos[v1] = vector3(0, 1, 0);
	aaPos[v2] = vector3(1, 0, 0);
	aaPos[v3] = vector3(-1, 0, 0);
	aaPos[v4] = vector3(0, 0.5, 1);
	//	create the tetrahedron
	Volume* v = *grid.create<Tetrahedron>(TetrahedronDescriptor(v1, v2, v3, v4));
	UG_LOG("cal tet" << endl);
	number vol = CalculateVolume(*v);
	BOOST_MESSAGE("volume of tetrahedron : " << vol);
}

BOOST_AUTO_TEST_SUITE_END();
