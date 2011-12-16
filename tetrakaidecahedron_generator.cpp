/*
 * tetrakaidekaeder_generator.cpp
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#include "tetrakaidecahedron_generator.h"
#include "generator.h"

#include "lib_grid/algorithms/extrusion/extrude.h"

using namespace ug;

namespace tkdGenerator {

/**
 * @param grid grid instance in which geometric objects will be created
 * @param positions
 * @param indices: numInds1, ind1_1, ind1_2, ..., numInds2, ind2_1, ind2_2, ...
 *
 * numInds == 4: tetrahedron
 * numInds == 5: pyramid
 * numInds == 6: prism
 * numInds == 8: hexahedron
 */
void createGrid(Grid& grid, const CoordsArray& positions,
		const IndexArray& indices) {

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
			aPosition, 10E-5);
}

void GenerateCorneocyteWithLipid(Grid& grid, number a_corneocyte,
		number w_corneocyte, number h_corneocyte, number d_lipid) {

//	// stolen parameters from old tkd modeler
//	number a1 = sqrt(
//			1.0 / 9.0 * h_corneocyte * h_corneocyte
//					+ 1.0 / 3.0 * pow((w_corneocyte - 2.0 * a_corneocyte), 2));
//
//	number alpha = acos((w_corneocyte - 2.0 * a_corneocyte) / (2.0 * a1));
//	number beta = 90.0 / 180.0 * PI
//			+ acos(1.0 / 3.0 * h_corneocyte / (a1 * sin(alpha)));
//	number gamma = acos(1.0 / 3.0 * h_corneocyte / a1) + 90.0 / 180.0 * PI;
//
//	number m1 = (d_lipid / 2) / tan(beta / 2);
//	number m2 = (d_lipid / 2) / tan(gamma / 2);
//
//	number a_lipid = (sqrt(3) + a_corneocyte + m1 + m2) / sqrt(3);
//	// fixme is this correct?
//	number h_lipid = h_corneocyte + d_lipid;
//	// fixme is this correct?
//	number w_lipid = w_corneocyte + d_lipid;
//
//	UG_LOG(
//			"a_c: " << a_corneocyte << "\n" << "w_c: " << w_corneocyte << "\n" << "h_c: " << h_corneocyte << endl);
//
//	UG_LOG(
//			"a_l: " << a_lipid << "\n" << "w_l: " << w_lipid << "\n" << "h_l: " << h_lipid << endl);

	CoordsArray posOut;
	IndexArray indsOut;

	Generator gCorneocyte(h_corneocyte, a_corneocyte, w_corneocyte, posOut,
			indsOut);

	gCorneocyte.createTKD();

	createGrid(grid, posOut, indsOut);

	std::vector<Face*> boundaryFaces;
	std::vector<VertexBase*> verticesOut;
	std::vector<EdgeBase*> edgesOut;
	for (FaceIterator iter = grid.begin<Face>();
		iter != grid.end<Face>(); ++iter)
	{
		if(IsVolumeBoundaryFace(grid, *iter))
			boundaryFaces.push_back(*iter);
	}
	UG_LOG("num faces: " << grid.num<Face>() << endl);
	UG_LOG("size boundaryFaces: " << boundaryFaces.size() << endl);

	Extrude(grid, &verticesOut, &edgesOut, &boundaryFaces, origin);

	//TODO scale verticesOut to match d_lipid
	// extruded tkd height = d_lipid + h_corneocyte
	number scale = 1;

	//	access vertex position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	for (size_t i = 0; i < verticesOut.size(); i++) {
//		VertexBase& vertex = *verticesOut[i];
//		v& coord = aaPos[i];
//		VecScale(coord, coord, scale);
	}
}

////////////////////////////////////////////////////////////////////////
///	test tetrakaidekahedron generator
void TestTKDGenerator(const char *outfile, number height, number baseEdgeLength,
		number diameter, number d_lipid) {
	Grid g(GRIDOPT_STANDARD_INTERCONNECTION);
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	g.attach_to_vertices(aPosition);

	GenerateCorneocyteWithLipid(g, baseEdgeLength, diameter, height, d_lipid);
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
