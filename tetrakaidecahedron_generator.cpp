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

		// all volumes belongs to corneocyte subset
		sh.assign_subset(vol, CORNEOCYTE);
	}

	// remove double vertices
	RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
			aPosition, 10E-5);
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
		aaPos[vrt] = vec;
	}

	// assign extruded selection to lipid subset
	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), LIPID);
	////////// End of extrusion

	// assigns wrong distanced faces to subset BROKEN
//	meassureLipidThickness(grid, sh, d_lipid, true);

	///// Correction of distance of side quadrilaterals
	// for all edges
	// fixme is it enough to iterate over subset LIPID only?
//	for (EdgeBaseIterator iter = sel.begin<EdgeBase>();
//			iter != sel.end<EdgeBase>(); iter++) {
//		EdgeBase* edge = *iter;
//		// for all boundary edges
//		if (IsBoundaryEdge3D(grid, edge)) {
//			vector<Face*> faces;
//			CollectAssociated(faces, grid, edge);
//			Face* triangles[2] = { NULL, NULL };
//			// for all edge associated faces
//			for (size_t i = 0; i < faces.size(); i++) {
//				Face* f = faces[i];
//				// assign faces if they are triangles
//				if (f->reference_object_id() == ROID_TRIANGLE) {
//					if (triangles[0])
//						triangles[1] = f;
//					else
//						triangles[0] = f;
//				}
//			}
//
//			// for all edges assiociated to 2 triangles
//			if (triangles[0] && triangles[1]) {
//				vector3 n1, n2;
//				CalculateNormal(n1, triangles[0], aaPos);
//				CalculateNormal(n2, triangles[1], aaPos);
//				number d = fabs(VecDot(n1, n2));
//				// ensure the 2 triangles does not have the same orientation
//				if (fabs(d - 1) > 10E-6) {
//					// faces[3] should be quadrilateral with wrong distance
//					UG_ASSERT(
//							faces[3]->reference_object_id() == ROID_QUADRILATERAL,
//							"faces[3] is not a quad");
//
//					/*for(sh.begin<Face>(BROKEN))
//					 *TODO ensure face[3] is broken...
//					 */
//
//					vector3 d1, d2;
//					vector3& v1 = aaPos[edge->vertex(0)], v2 =
//							aaPos[edge->vertex(1)];
//
//					VecSubtract(d1, v1, v2);
//					VecScale(d2, d1, -1);
//					vector3 offset;
//					vector3 normal;
//					CalculateNormal(normal, faces[3], aaPos);
//					number len1 = VecDot(normal, d1);
//					number len2 = VecDot(normal, d2);
//					// direction * normal = 0
//					UG_LOG("n * d1: " << len1 << "\t n*d2: " << len2 << endl);
//					// TODO move v1, v2 in direction of center of edge
////					VecScale(v1, v1, 1. / len2);
////					VecScale(v2, v2, 1. / len1);
//					sh.assign_subset(edge, 3);
//				}
//			}
//		}
//	}

	Selector selEdges(grid);

	vector<Face*> faces;
	for (EdgeBaseIterator iter = grid.begin<Edge>(); iter != grid.end<Edge>();
			iter++) {
		EdgeBase* edge = *iter;
		if (IsBoundaryEdge3D(grid, edge)) {
			CollectAssociated(faces, grid, edge, true);

			Face* f[2] = { NULL, NULL };
			for (uint i = 0; i < faces.size(); i++) {
				if (LiesOnBoundary(grid, faces[i])) {
					if (f[0])
						f[1] = faces[i];
					else
						f[0] = faces[i];
				}
			}

			vector3 n1, n2;
			CalculateNormal(n1, faces[0], aaPos);
			CalculateNormal(n2, faces[1], aaPos);

			number d1 = fabs(VecDot(n1, n2));
			number epsilon = 10E-6;
			number delta = fabs(d1 - 1);
			UG_LOG("d1: " << d1 << "\tdelta: " << delta << endl);
			if (delta > epsilon) {
//				UG_LOG("selecting edge..." << endl);
				selEdges.select(edge);
			} else {
				UG_LOG("normals are parallel." << endl);
			}
		}
	}

	vector<EdgeBase*> face_egdes;
	for (FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>();
			iter++) {
		Face* face = *iter;
		if (IsBoundaryFace3D(grid, face)) {
			CollectAssociated(face_egdes, grid, face, true);

			int selected = 0, num_edges = face->num_edges();
			for (int i = 0; i < num_edges; i++) {
				if (selEdges.is_selected(grid.get_edge(face, i)))
					selected++;
			}

			if (selected == num_edges) {
				selEdges.deselect(face_egdes.begin(), face_egdes.end());
			}
		}
	}

	sh.assign_subset(selEdges.begin<Edge>(), selEdges.end<Edge>(), 3);

	sh.subset_info(3).name = "scale_edges";
	sh.subset_info(3).color = vector4(1, 0, 0, 0);

	//// Copy extruded tkd in each dimension
	uint count = rows * cols * high;
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
		// grid.erase(sel.begin<Volume>(), sel.end<Volume>());

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
