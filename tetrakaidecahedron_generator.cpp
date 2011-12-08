/*
 * tetrakaidekaeder_generator.cpp
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#include "tetrakaidecahedron_generator.h"
#include "generator.h"

using namespace ug;

namespace tkdGenerator {

void GenerateTetrakaidecahedron(Grid& grid, number& height,
		number& baseEdgeLength, number& diameter) {

	CoordsArray positions;
	IndexArray indices;
	Generator generator(height, baseEdgeLength, diameter, positions, indices);

//	fill the arrays
	generator.createTKD();

//	access vertex position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	generate vertices in the grid and store them in an array so that we can index them
	std::vector<VertexBase*> vertices(positions.size());

	for (size_t i = 0; i < positions.size(); ++i) {
		VertexBase* v = *grid.create<Vertex>();
		aaPos[v] = positions[i];
		vertices[i] = v;
	}

//	the VolumeDescriptor will be used to create new volumes
	VolumeDescriptor vd;

//	create the elements from the given indices
	for (size_t i = 0; i < indices.size();) {
		int num = indices[i++];
		vd.set_num_vertices(num);
		for (int j = 0; j < num; ++j)
			vd.set_vertex(j, vertices[indices[i++]]);

		switch (num) {
		case 4:
			grid.create<ug::Tetrahedron>(vd);
			break;
		case 5:
			grid.create<ug::Pyramid>(vd);
			break;
		case 6:
			grid.create<ug::Prism>(vd);
			break;
		case 8:
			grid.create<ug::Hexahedron>(vd);
			break;
		}
	}

	// remove double vertices
	RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
			aPosition, 0.1);
}

void GenerateCorneocyteWithLipid(number a_corneocyte, number width, number H,
		number d_lipid) {
// stolen parameters from old tkd modeler

//	number a1 = sqrt(
//			1.0 / 9.0 * H * H
//					+ 1.0 / 3.0 * pow((width - 2.0 * a_corneocyte), 2));
//
//	number alpha = acos((width - 2.0 * a_corneocyte) / (2.0 * a1));
//	number beta = 90.0 / 180.0 * PI + acos(1.0 / 3.0 * H / (a1 * sin(alpha)));
//	number gamma = acos(1.0 / 3.0 * H / a1) + 90.0 / 180.0 * PI;
//
//	number m1 = (d_lipid / 2) / tan(beta / 2);
//	number m2 = (d_lipid / 2) / tan(gamma / 2);
//
//	number a_lipid = (sqrt(3) + a_corneocyte + m1 + m2) / sqrt(3);
//	number h_lipid = d_lipid + H;
//	number w_lipid;

	// create corneocyte and lipid tkd and merge them somehow
	//GenerateTetrakaidecahedron()
}

////////////////////////////////////////////////////////////////////////
///	test tetrakaidekahedron generator
void TestTKDGenerator(const char *outfile, number height, number baseEdgeLength,
		number diameter) {
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);
	tkdGenerator::GenerateTetrakaidecahedron(g, height, baseEdgeLength,
			diameter);
	SaveGridToFile(g, sh, outfile);
}

extern "C" void InitUGPlugin(ug::bridge::Registry* reg,
		std::string parentGroup) {
	std::string grp(parentGroup);
	grp.append("tkd_generator/");

	//	add TKD-Generator method
	reg->add_function("TestTKDGenerator", &TestTKDGenerator, grp);
}

} // end of namespace tkdGenerator
