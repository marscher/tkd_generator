/*
 * tetrakaidekaeder_generator.cpp
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#include "tetrakaidecahedron_generator.h"
#include "generator.h"

#include "lib_grid/lib_grid.h"
#include "registry/registry.h"

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
void createGrid(Grid& grid, SubsetHandler& sh, const CoordsArray& positions,
		const IndexArray& indices) {

	// access vertex position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	// generate vertices in the grid and store them in an array so that we can index them
	std::vector<VertexBase*> vertices(positions.size());

	for (size_t i = 0; i < positions.size(); ++i) {
		VertexBase* v = *grid.create<Vertex>();
		aaPos[v] = positions[i];
		vertices[i] = v;
	}

	// the VolumeDescriptor will be used to create new volumes
	VolumeDescriptor vd;

	// create the elements from the given indices
	for (size_t i = 0, count = 0; i < indices.size(); count++) {
		int num = indices[i++];
		vd.set_num_vertices(num);
		for (int j = 0; j < num; ++j)
			vd.set_vertex(j, vertices[indices[i++]]);

		Volume* vol;
		switch (num) {
		case 4:
			vol = *grid.create<ug::Tetrahedron>(vd);
			break;
		case 5:
			vol = *grid.create<ug::Pyramid>(vd);
			break;
		case 6:
			vol = *grid.create<ug::Prism>(vd);
			break;
		case 8:
			vol = *grid.create<ug::Hexahedron>(vd);
			break;
		}

		// first 63 volumes belongs to corneocyte subset, rest to lipid
		if (count < 63) {
			sh.assign_subset(vol, CORNEOCYTE);
		} else {
			sh.assign_subset(vol, LIPID);
		}
	}

	// remove double vertices
	RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
			aPosition, 10E-5);
}

void GenerateCorneocyteWithLipid(Grid& grid, SubsetHandler& sh,
		number a_corneocyte, number w_corneocyte, number h_corneocyte,
		number d_lipid, int rows, int cols, int high) {

	grid.set_options(GRIDOPT_STANDARD_INTERCONNECTION);
	sh.assign_grid(grid);
	// Subset 1 for corneocytes (same as in Feuchters tkdmodeller)
	sh.set_default_subset_index(CORNEOCYTE);
	sh.subset_info(CORNEOCYTE).name = "corneocytes";
	sh.subset_info(LIPID).name = "lipid matrix";
	sh.subset_info(BROKEN).name = "broken";

	// argb: green
	sh.subset_info(LIPID).color = vector4(0, 1, 0, 0);
	// argb: blue
	sh.subset_info(CORNEOCYTE).color = vector4(0, 0, 1, 0);

	grid.attach_to_vertices(aPosition);

	TKDGeometryGenerator generator(h_corneocyte, a_corneocyte, w_corneocyte, d_lipid);
	generator.createDomain();
	createGrid(grid, sh, generator.getPositions(), generator.getIndices());

	int count = rows * cols * high;

	// copy first tkd count times with proper offsets
	if (count > 1) {
		UG_LOG("creating " << count << " cells with lipid matrix." << endl);

		Selector sel(grid);
		sel.enable_autoselection(false);
		// select lipid and corneocyte of first tkd for duplication
		sel.select(grid.begin<Volume>(), grid.end<Volume>());

		// do not select duplicates so we can continue to work with initial selection
		bool selectNew = false;
		bool deselectOld = false;

		number offsetX = 0;
		number offsetY = 0;
		number offsetZ = 0;

		// 1 tkd = three segments with each 1/3 h
		number third_height = (1. / 3 * generator.getHeight());
		number s_lipid = generator.getOverlap();
		// overlap in x and y direction is overlap of lipid
		for (int i = 0; i < rows; i++, offsetX += s_lipid) {
			offsetY = 0;
			// every even row is shifted by 1/3 h_lipid
			if (i % 2) {
				offsetZ += third_height;
			}
			for (int j = 0; j < cols; j++, offsetY += s_lipid) {
				offsetZ = 0;
				// every even column is shifted by 1/3 h_lipid
				if (j % 2) {
					offsetZ += third_height;
				}
				for (int k = 0; k < high;
						k++, offsetZ += generator.getHeight()) {
					vector3 offset(offsetX, offsetY, offsetZ);
					Duplicate(grid, sel, offset, aPosition, deselectOld,
							selectNew);
				}
			}
		}

// todo delete first created tkd, because it is offset o the others
//		grid.erase(sel.begin<Volume>(), sel.end<Volume>());

// finally remove double vertices
//		RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
//				aPosition, 10E-5);
	}
}

////////////////////////////////////////////////////////////////////////
///	test tetrakaidekahedron generator
void TestTKDGenerator(const char *outfile, number height, number baseEdgeLength,
		number diameter, number d_lipid, int rows, int cols, int high) {

	Grid g;
	SubsetHandler sh;

	GenerateCorneocyteWithLipid(g, sh, baseEdgeLength, diameter, height,
			d_lipid, rows, cols, high);

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
