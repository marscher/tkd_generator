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

void GenerateCorneocyteWithLipid(Grid& grid, SubsetHandler& sh,
		number a_corneocyte, number w_corneocyte, number h_corneocyte,
		number d_lipid, int rows, int cols, int high) {

	// stolen parameters from old tkd modeler
	number a1 = sqrt(
			1.0 / 9.0 * h_corneocyte * h_corneocyte
					+ 1.0 / 3.0 * pow((w_corneocyte - 2.0 * a_corneocyte), 2));

	number alpha = acos((w_corneocyte - 2.0 * a_corneocyte) / (2.0 * a1));
	number beta = 90.0 / 180.0 * M_PI
			+ acos(1.0 / 3.0 * h_corneocyte / (a1 * sin(alpha)));
	number gamma = acos(1.0 / 3.0 * h_corneocyte / a1) + 90.0 / 180.0 * M_PI;

	number m1 = (d_lipid / 2) / tan(beta / 2);
	number m2 = (d_lipid / 2) / tan(gamma / 2);

	number a_lipid = (sqrt(3) + a_corneocyte + m1 + m2) / sqrt(3);
	///////// end of stolen stuff

	// extruded tkd height = d_lipid + h_corneocyte
	number h_lipid = h_corneocyte + d_lipid;
	number h_scale = h_lipid / h_corneocyte;
	UG_LOG("h_scale: " << h_scale << endl);

	// fixme this should be the right value!
	number w_lipid = w_corneocyte + d_lipid;

	// overlap of outer lipid tkd
	number s_lipid = 1 / sqrt(3) * (a_lipid - 2 * w_lipid);

	UG_LOG(
			"a_c: " << a_corneocyte << "\n" << "w_c: " << w_corneocyte << "\n" << "h_c: " << h_corneocyte << endl);

	UG_LOG(
			"a_l: " << a_lipid << "\n" << "w_l: " << w_lipid << "\n" << "h_l: " << h_lipid << endl);

	Generator gCorneocyte(h_corneocyte, a_corneocyte, w_corneocyte);
	gCorneocyte.createTKD();
	createGrid(grid, gCorneocyte.getPosOut(), gCorneocyte.getIndsOut());

	std::vector<Face*> boundaryFaces;
	std::vector<Face*>::iterator bfIter;

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

	// access distances from inner faces to outer faces
	ANumber aDist;
	grid.attach_to_faces(aDist);
	Grid::FaceAttachmentAccessor<ANumber> aaDist(grid, aDist);

	// scale every vertex of selection
	for (VertexBaseIterator iter = sel.begin<VertexBase>();
			iter != sel.end<VertexBase>(); iter++) {
		VertexBase* vrt = *iter;
		vector3& vec = aaPos[vrt];
		VecScale(vec, vec, h_scale);
		aaPos[vrt] = vec;
	}

	// assign to lipid subset
	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), 1);

	// log element height
	sel.clear();
	std::vector<Face*> faces;
	std::vector<Face*> toCorrect;
	SelectAreaBoundaryFaces(sel, sh.begin<Volume>(1), sh.end<Volume>(1));

	for (VolumeIterator iter = sh.begin<Volume>(1); iter != sh.end<Volume>(1);
			iter++) {
		Volume* v = *iter;
		CollectAssociated(faces, grid, v);
		Face* selFaces[2] = { NULL, NULL };
		for (size_t i_face = 0; i_face < faces.size(); ++i_face) {
			Face* f = faces[i_face];
			if (sel.is_selected(f)) {
				if (selFaces[0]) {
					selFaces[1] = f;
				} else {
					selFaces[0] = f;
				}
			}
		}

		UG_ASSERT(selFaces[0] && selFaces[1],
				"There should be exactly 2 selected faces!");

		vector3 center1 = CalculateCenter(selFaces[0], aaPos);
		vector3 center2 = CalculateCenter(selFaces[1], aaPos);
		vector3 normal;
		CalculateNormal(normal, selFaces[1], aaPos);

		vector3 point;
		ProjectPointToPlane(point, center1, center2, normal);
		number dist = VecDistance(center1, point);
		UG_LOG("distance : " << dist << endl);
		number diff = dist - d_lipid / 2;
		faces.clear();
		if (fabs(diff > 0.1)) {
			sel.select(v);
			CollectAssociated(faces, grid, v);
			for (size_t i_face = 0; i_face < faces.size(); ++i_face) {
				Face* f = faces[i_face];
				if (IsVolumeBoundaryFace(grid, f)) {
					aaDist[f] = diff;
					toCorrect.push_back(f);
					UG_LOG("aadist[f]: " << aaDist[f] << endl);
				}
			}
		}
	}

//	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), 2);
	UG_LOG("beginning correction..." << endl);
	for (size_t i = 0; i < toCorrect.size(); i++) {
		Face* outerFace = toCorrect[i];
		number correction = aaDist[outerFace];
		for (size_t j = 0; j < outerFace->num_vertices(); j++) {
			VertexBase* vrt = outerFace->vertex(j);
			vector3& vec = aaPos[vrt];
			//fixme
//			vec -= correction;
			aaPos[vrt] = vec;
		}
	}

	// detach distance information from faces because is not needed anymore
	grid.detach_from_faces(aDist);

	int count = rows * cols * high;
	UG_LOG("creating " << count << " tkds." << endl);

	// copy first tkd count times with proper offsets
	if (count > 1) {
		// clear current selection
		sel.clear();
		// select lipid and corneocyte of first tkd for duplication
		sel.select(grid.begin<Volume>(), grid.end<Volume>());

		// do not select duplicates so we can continue to work with initial selection
		bool selectNew = false;
		bool deselectOld = false;

		number offsetX = 0;
		number offsetY = 0;
		number offsetZ = 0;

		// 1 tkd = three segments with each 1/3 h
		const number third_height = (1. / 3 * h_lipid);

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
				for (int k = 0; k < high; k++, offsetZ += h_lipid) {
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
	Grid g(GRIDOPT_STANDARD_INTERCONNECTION);
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);
	sh.subset_info(0).name = "corneocytes";
	sh.subset_info(1).name = "lipid matrix";

	// argb: green
	sh.subset_info(0).color = vector4(0, 1, 0, 0);
	// argb: red
	sh.subset_info(1).color = vector4(0, 0, 1, 0);

	g.attach_to_vertices(aPosition);

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
