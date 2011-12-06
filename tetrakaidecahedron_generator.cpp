/*
 * tetrakaidekaeder_generator.cpp
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#include "tetrakaidecahedron_generator.h"
#include "geometric_helper.h"

using namespace ug;

namespace tkdGenerator {

void GenerateTetrakaidecahedron(Grid& grid, number& height,
		number& baseEdgeLength, number& diameter) {

	CoordsArray positions;
	IndexArray indices;

//	fill the arrays
	GenerateTetrakaidecahedron(positions, indices, height, baseEdgeLength,
			diameter);

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
			grid.create<Tetrahedron>(vd);
			break;
		case 5:
			grid.create<Pyramid>(vd);
			break;
		case 6:
			grid.create<Prism>(vd);
			break;
		case 8:
			grid.create<Hexahedron>(vd);
			break;
		}
	}

	// remove double vertices
	RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
			aPosition, 0.1);
}

void GenerateCorneocyteWithLipid(number a_corneocyte, number width, number H,
		number d_lipid) {

	number a1 = sqrt(
			1.0 / 9.0 * H * H
					+ 1.0 / 3.0 * pow((width - 2.0 * a_corneocyte), 2));

	number alpha = acos((width - 2.0 * a_corneocyte) / (2.0 * a1));
	number beta = 90.0 / 180.0 * PI + acos(1.0 / 3.0 * H / (a1 * sin(alpha)));
	number gamma = acos(1.0 / 3.0 * H / a1) + 90.0 / 180.0 * PI;

	number m1 = (d_lipid / 2) / tan(beta / 2);
	number m2 = (d_lipid / 2) / tan(gamma / 2);

	number a_lipid = (sqrt(3) + a_corneocyte + m1 + m2) / sqrt(3);
	number h_lipid = d_lipid + H;
	number w_lipid;

	// create corneocyte and lipid tkd and merge them somehow
	//GenerateTetrakaidecahedron()
}

/**
 *
 */
void GenerateTetrakaidecahedron(CoordsArray& posOut, IndexArray& indsOut,
		number& height, number& baseEdgeLength, number& diameter) {
	// baseEdgeLength
	number a = baseEdgeLength;
	// height of one level of decomposition which is 1/3 of overall height of the tkd
	number h = height / 3;
	// quantity s, overlap of two aligned tkd's
	number s = 1 / sqrt(3) * (diameter - 2 * baseEdgeLength);
	// height of base triangle of top inner prism
	number g = sqrt(3) * a / 2;
	// height of base triangle of tetrahedron of ObenAussenPr2T
	number b = sqrt(3) * s / 2;

	UG_LOG("g: " << g << endl);
	UG_LOG("s: " << s << endl);
	UG_LOG("b: " << b << endl);

	// origin of construction of tkd
	const v origin(0, 0, 0); //(10, 10, 10);

	// Rotation matrix with initial rotation angle of 0 (rad)
	rotationMatrix R(0);

	// create G(Ki -> ObenInnen) = 6 prism with equilateral sites
	// rotate prism around (0, y, 0) with step angle 60° = PI/3 (rad)
	for (number t = 0; t < 2 * PI; t += PI / 3) {
		// switch angle of rotation matrix
		R.setAngle(t);

		// begin in origin
		v v1 = myTransform(v(0, 0, 0), R, origin);
		// create surrounding vertices relative to origin
		v v2 = myTransform(v(a, 0, 0), R, origin);
		v v3 = myTransform(v(a / 2, 0, g), R, origin);
		v v4 = myTransform(v(0, h, 0), R, origin);
		v v5 = myTransform(v(a, h, 0), R, origin);
		v v6 = myTransform(v(a / 2, h, g), R, origin);

		createPrism(v1, v2, v3, v4, v5, v6, posOut, indsOut);
	}

	// create G(Ki -> ObenAussenPr2T)
	for (number t = 2. / 3 * PI; t < 2 * PI; t += (2. / 3 * PI)) {
		R.setAngle(t);
		 b = sqrt(2) * s;
		// TODO b ist nicht das korrekte offset, da es nicht paralell zur x achse liegt.
		// left tetrahedron right of prism of ObenAussenPr
		v v1_a = myTransform(v(a / 2, 0, g), R, origin);
		v v2_a = myTransform(v(a / 2, 0, g + s), R, origin);
		v v3_a = myTransform(v(a / 2 + b, 0, g + s / 2), R, origin);
		// top
		v v4_a = myTransform(v(a / 2, h, g), R, origin);
		createTetrahedron(v1_a, v2_a, v3_a, v4_a, posOut, indsOut);

//		// Prism of ObenAussenPr2T
		//TODO tiefe von prisma ist s/2 wie von tetraeder, damit die auf einer ebene liegen
		v v1p = myTransform(v(a / 2, 0, g), R, origin);
		v v2p = myTransform(v(a / 2, h, g), R, origin);
		v v3p = myTransform(v(a / 2 + b, 0, g + s / 2), R, origin);

		v v4p = myTransform(v(a, 0, 0), R, origin);
		v v5p = myTransform(v(a, h, 0), R, origin);
		v v6p = myTransform(v(a + b, 0, s / 2), R, origin);

		createPrism(v1p, v2p, v3p, v4p, v5p, v6p, posOut, indsOut);

		// right of prism ObenAussenPr2T
		v v1_b = myTransform(v(a, 0, 0), R, origin);
		v v2_b = myTransform(v(a, b, 0), R, origin);
		v v3_b = myTransform(v(a + b, 0, 0), R, origin);

		v v4_b = myTransform(v(a, h, 0), R, origin);
//		createTetrahedron(v1_b, v2_b, v3_b, v4_b, posOut, indsOut), R, origin);
//
		// create G(Ki -> ObenAussenPr):
		// create 3 prisms every 120° or 2/3 PI (rad)
		v v1 = myTransform(v(a / 2, 0, g), R, origin);
		v v2 = myTransform(v(a / 2, h, g), R, origin);
		v v3 = myTransform(v(a / 2, 0, g + s), R, origin);
		v v4 = myTransform(v(-a / 2, 0, g), R, origin);
		v v5 = myTransform(v(-a / 2, h, g), R, origin);
		v v6 = myTransform(v(-a / 2, 0, g + s), R, origin);

		createPrism(v1, v2, v3, v4, v5, v6, posOut, indsOut);
	}
}

////////////////////////////////////////////////////////////////////////
///	test tetrakaidekahedron generator
void TestTKDGenerator(const char* outfile, number height, number baseEdgeLength,
		number diameter) {
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);

	g.attach_to_vertices(aPosition);

	tkdGenerator::GenerateTetrakaidecahedron(g, height, baseEdgeLength,
			diameter);
	SaveGridToFile(g, sh, outfile);
}

extern "C" UG_API void InitUGPlugin(ug::bridge::Registry* reg,
		std::string parentGroup) {
	std::string grp(parentGroup);
	grp.append("tkd_generator/");

	//	add TKD-Generator method
	reg->add_function("TestTKDGenerator", &TestTKDGenerator, grp);
}

} // end of namespace tkdGenerator
