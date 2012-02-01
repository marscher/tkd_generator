/*
 * tetrakaidekaeder_generator.cpp
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#include "tetrakaidecahedron_generator.h"
#include "generator.h"
#include "test/testHelper.h"

#include "lib_grid/lib_grid.h"
#include "registry/registry.h"

using namespace ug;

namespace tkdGenerator {

/**
 * creates lib_grid objects in given grid reference according to given positions and indices
 * @param grid grid instance in which geometric objects will be created
 * @param positions
 * @param indices: numInds1, ind1_1, ind1_2, ..., numInds2, ind2_1, ind2_2, ...
 *
 * numInds == 4: tetrahedron
 * numInds == 5: pyramid
 * numInds == 6: prism
 * numInds == 8: hexahedron
 */
void createGridFromArrays(Grid& grid, SubsetHandler& sh,
		const CoordsArray& positions, const IndexArray& indices) {

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

//	sh.assign_subset(grid.begin<Vertex>(), grid.end<Vertex>(), CORNEOCYTE);
//	sh.assign_subset(grid.begin<Edge>(), grid.end<Edge>(), CORNEOCYTE);
	sh.assign_subset(grid.begin<Volume>(), grid.end<Volume>(), CORNEOCYTE);
}

void GenerateCorneocyteWithLipid(Grid& grid, SubsetHandler& sh,
		number a_corneocyte, number w_corneocyte, number h_corneocyte,
		number d_lipid, int rows, int cols, int high) {

	//// set grid options
	grid.set_options(GRIDOPT_STANDARD_INTERCONNECTION);
	grid.attach_to_vertices(aPosition);

	//// Subset informations
	sh.assign_grid(grid);
	// Subset 1 for corneocytes (same as in Feuchters tkdmodeller)
	sh.set_default_subset_index(CORNEOCYTE);
	sh.subset_info(CORNEOCYTE).name = "corneocytes";
	sh.subset_info(LIPID).name = "lipid matrix";
	sh.subset_info(BROKEN).name = "broken";
	//// Colors
	// argb: green
	sh.subset_info(LIPID).color = vector4(0, 1, 0, 0);
	// argb: blue
	sh.subset_info(CORNEOCYTE).color = vector4(0, 0, 1, 0);

	//// create coordinates and vertex indices for 1 tetrakaidecahedron
	TKDGeometryGenerator generator(h_corneocyte, a_corneocyte, w_corneocyte,
			d_lipid);
	generator.createDomain();
	//// fill the grid object
	createGridFromArrays(grid, sh, generator.getPositions(),
			generator.getIndices());

	//// extrude the tkd to match given lipid thickness size and correct quadrilaterals afterwards
	// extruded tkd height = d_lipid + h_corneocyte
	number h_lipid = h_corneocyte + d_lipid;
	number h_scale = h_lipid / h_corneocyte;
//	UG_LOG("h_scale: " << h_scale << endl);

	vector<Face*> boundaryFaces;
	for (FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>();
			++iter) {
		if (IsVolumeBoundaryFace(grid, *iter))
			boundaryFaces.push_back(*iter);
	}

	// select vertices which are beeing extruded
	Selector sel(grid);
	sel.enable_autoselection(true);

	// extruding boundary faces
	Extrude(grid, NULL, NULL, &boundaryFaces, origin);

	//	access vertex position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	// scale every vertex of selection
	for (VertexBaseIterator iter = sel.begin<VertexBase>();
			iter != sel.end<VertexBase>(); iter++) {
		VertexBase* vrt = *iter;
		vector3& vec = aaPos[vrt];
		VecScale(vec, vec, h_scale);
	}

	// assign extruded selection to lipid subset
	// todo edges for testing only

	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), LIPID);
	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), LIPID);
//	for (EdgeIterator iter = sel.begin<Edge>(); iter != sel.end<Edge>();
//			iter++) {
//		EdgeBase* e  = *iter;
//		if (LiesOnBoundary(grid, e)) {
//			sh.assign_subset(e, LIPID);
//		}
//	}
	sel.clear();
	////////// End of extrusion

	// assigns wrong distanced faces to subset BROKEN
//	meassureLipidThickness(grid, sh, d_lipid, true);

	///// Correction of distance of side quadrilaterals

//	vector<Face*> faces;
//	for (EdgeBaseIterator iter = grid.begin<Edge>(); iter != grid.end<Edge>();
//			iter++) {
//		EdgeBase* edge = *iter;
//		if (IsBoundaryEdge3D(grid, edge)) {
//			CollectAssociated(faces, grid, edge, true);
//
//			Face* f[2] = { NULL, NULL };
//			for (uint i = 0; i < faces.size(); i++) {
//				if (LiesOnBoundary(grid, faces[i])) {
//					if (f[0])
//						f[1] = faces[i];
//					else
//						f[0] = faces[i];
//				}
//			}
//
//			UG_ASSERT(f[0] && f[1], "two faces are needed.");
//
//			vector3 n1, n2;
//			CalculateNormal(n1, f[0], aaPos);
//			CalculateNormal(n2, f[1], aaPos);
//
//			number d1 = fabs(VecDot(n1, n2));
//			if (d1 > .5) {
//				sel.select(edge);
//			}
//		}
//	}
//
//	vector<EdgeBase*> face_egdes;
//	for (FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>();
//			iter++) {
//		Face* face = *iter;
//		if (IsBoundaryFace3D(grid, face)) {
//			CollectAssociated(face_egdes, grid, face, true);
//
//			int selected = 0, num_edges = face->num_edges();
//			for (int i = 0; i < num_edges; i++) {
//				if (sel.is_selected(grid.get_edge(face, i)))
//					selected++;
//			}
//
//			if (selected == num_edges) {
//				sel.deselect(face_egdes.begin(), face_egdes.end());
//			}
//		}
//	}
//
//	UG_LOG("selected edges: " << sel.num<Edge>() << endl);
//	sh.assign_subset(sel.begin<Edge>(), sel.end<Edge>(), 3);
//
//	sh.subset_info(3).name = "scale_edges";
//	sh.subset_info(3).color = vector4(1, 0, 0, 0);

	//// Copy extruded tkd in each dimension
	uint count = rows * cols * high;
	if (count > 1) {
		UG_LOG("creating " << count << " cells with lipid matrix." << endl);

		sel.clear();
		sel.enable_autoselection(false);
		// select lipid and corneocyte of first tkd for duplication
		sel.select(grid.begin<Volume>(), grid.end<Volume>());

		// do not select duplicates so we can continue to work with initial selection
		bool selectNew = false;
		bool deselectOld = false;

		vector3 offset(0, 0, 0);

		// determine dimensions of extruded tkd
		vector3 min, max, dimension;
		CalculateBoundingBox(min, max, grid.begin<Vertex>(), grid.end<Vertex>(),
				aaPos);

		VecSubtract(dimension, max, min);
		UG_LOG(
				"bbox: " << min << " "<< max << "\ndimensions: " << dimension << endl);

		// fixme x != y
		number d = dimension.x;
		number h = dimension.z;
		// determine length of 1 edge of base hexagon (equilateral triangles)
		number a = a_corneocyte * h_scale;

//		number w = 1. / 6
//				* (sqrt(3) * sqrt(-9 * a * a + 9 * d * d - h * h) - 3 * a);
		number w = w_corneocyte * h_scale;//(sqrt(3) * a + 3 * dimension.x) / 2 * sqrt(3);
		number s = 1 / sqrt(3) * (w - 2 * a);
		UG_LOG("a: "<< a << " w : " << w << " s: " << s << endl);
		UG_ASSERT(w > 2*a, "geometrical constraints not met.");

		// 1 tkd = three segments with each 1/3 h

		number h_shift = 2 * h / 3;
		for (int i = 0; i < rows; i++) {
			offset.x += 2*s;
			// reset y offset, as we are in a new row
			offset.y = 0;
			// every even row is shifted by 2/3 h_lipid
			if (i % 2) {
				offset.z += h_shift;
			}

			for (int j = 0; j < cols; j++) {
				offset.y += d;
				// reset z offset, as we are in a new col
				offset.z = 0;
				// every even column is shifted by 2/3 h_lipid
				if (j % 2) {
					offset.z += h_shift;
				}

				for (int k = 0; k < high; k++) {
					offset.z += h;
					UG_LOG(
							"offset(" << i << ", " << j << ", " << k << "): " << offset << endl);
					Duplicate(grid, sel, offset, aPosition, deselectOld,
							selectNew);
				}
			}
		}

		// delete first created tkd, because it is offset o the others
		grid.erase(sel.begin<Vertex>(), sel.end<Vertex>());
		grid.erase(sel.begin<Edge>(), sel.end<Edge>());
		grid.erase(sel.begin<Volume>(), sel.end<Volume>());

		// finally remove double vertices
		RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
				aPosition, 10E-5);
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

/**
 * to call from lua...
 */
void createTKDDomain(Grid& grid, number height, number baseEdgeLength,
		number diameter, number d_lipid, int rows, int cols, int high) {

	SubsetHandler sh;
	GenerateCorneocyteWithLipid(grid, sh, baseEdgeLength, diameter, height,
			d_lipid, rows, cols, high);
}

extern "C" void InitUGPlugin(ug::bridge::Registry* reg,
		std::string parentGroup) {
	std::string grp(parentGroup);
	grp.append("tkd_generator/");

	//	add TKD-Generator method
	reg->add_function("TestTKDGenerator", &TestTKDGenerator, grp);
	reg->add_function("createTKDDomain", &createTKDDomain, grp);
}

} // end of namespace tkdGenerator
