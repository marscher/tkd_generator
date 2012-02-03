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

#include <set>
//#include <map>

using namespace ug;
using namespace std;

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

	sh.assign_subset(grid.begin<Volume>(), grid.end<Volume>(), CORNEOCYTE);
}

struct vecComperator {
	// returns a < b
	bool operator()(const vector3& a, const vector3& b) const {
		number SMALL = 10E-6;
		if (a.x < b.x - SMALL)
			return true;
		if (a.x > b.x + SMALL)
			return false;

		if (a.y < b.y - SMALL)
			return true;
		if (a.y > b.y + SMALL)
			return false;

		if (a.z < b.z - SMALL)
			return true;
		if (a.z > b.z + SMALL)
			return false;

		return false;
	}
};

vector<pair<VertexBase*, vector3> > generateLipidMatrixForSingleTKD(Grid& grid,
		SubsetHandler& sh, number d_lipid) {

	vector<Face*> faces;
	set<vector3, vecComperator> normals;
	set<vector3, vecComperator>::iterator normalsIter;

	// stores pairs of vertices with their final shifts
	vector<pair<VertexBase*, vector3> > vertexShifts;

	// select boundaryFaces for extrusion
	for (FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>();
			++iter) {
		if (IsVolumeBoundaryFace(grid, *iter))
			faces.push_back(*iter);
	}

	// auto select vertices which are beeing extruded
	Selector sel(grid);
	sel.enable_autoselection(true);

	// extruding boundary faces
	Extrude(grid, NULL, NULL, &faces, origin);

	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), LIPID);

	//	access vertex position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	// for all extruded vertices
	for (VertexBaseIterator iter = sel.begin<VertexBase>();
			iter != sel.end<VertexBase>(); iter++) {

		VertexBase* v = *iter;
		// for all boundary vertices
		if (IsBoundaryVertex3D(grid, v)) {
			CollectAssociated(faces, grid, v);
			// collect normals only for this vertex
			normals.clear();
			// for all faces associated to v
			for (size_t i = 0; i < faces.size(); i++) {
				if (IsBoundaryFace3D(grid, faces[i])) {
					vector3 normal;
					CalculateNormal(normal, faces[i], aaPos);
					normals.insert(normal);
				}
			}

			uint num_normals = normals.size();
			normalsIter = normals.begin();

			switch (num_normals) {
			// vertex lies in bound of a tkd side hexagon
			case 1: {
				vector3 n = *normalsIter;
				// scale normal so it fits to half of lipid thickness
				VecScale(n, n, d_lipid / 2);
				// make copy of vertex position attachment
				vector3 a = aaPos[v];
				// shift a by scaled normal of face
				VecAdd(a, a, n);
				// store shift vector p to corresponding vertex v
				vertexShifts.push_back(make_pair(v, a));
				break;
			}

			// vertex lies on an edge between an hexagon and a quadrilateral
			case 3: {
				sh.assign_subset(v, 2);

				// init normals of vertex v associated faces
				vector3 na = *normalsIter++;
				vector3 nb = *normalsIter++;
				vector3 nz = *normalsIter;

				// a: vertex position
				// b: support vector of plane nb
				// z: support vector of plane nz
				// nc: lies in plane a
				// nt: lies in plane b?
				// c: intersection of nt with plane b
				// p: intersection of of nc with plane z
				vector3 a = aaPos[v], b, z, nc, nt, c, p;

				// init support vectors b, z
				VecAdd(b, a, nb);
				VecAdd(z, a, nz);

//				UG_LOG("na: " << na << "\tnb: "<< nb << "\tnz: "<< nz << endl);

				// scale normals to length d_lipid/2
				VecScale(na, na, d_lipid / 2);
				VecScale(nb, nb, d_lipid / 2);
				VecScale(na, na, d_lipid / 2);

				// calculate vectors lying in planes na, nc
				VecCross(nc, na, nb);
				VecCross(nt, na, nc);

				number tmp;
				RayPlaneIntersection(c, tmp, a, nt, b, nb);
				RayPlaneIntersection(p, tmp, c, nc, z, nz);

				// store shift vector p to corresponding vertex v
				vertexShifts.push_back(make_pair(v, p));
				break;
			}
			default:
				UG_ASSERT(true, "should never get here!");
			}
		} // end of switch
	} // end of calculate shifts for

	// apply shifts to all vertices
	for (size_t i = 0; i < vertexShifts.size(); i++) {
		VertexBase* v = vertexShifts[i].first;
		vector3& shift = vertexShifts[i].second;
		vector3& pos = aaPos[v];

		VecAdd(pos, pos, shift);
	}

	return vertexShifts;
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

	vector<pair<VertexBase*, vector3> > shifts =
			generateLipidMatrixForSingleTKD(grid, sh, d_lipid);

	// assigns wrong distanced faces to subset BROKEN
//	meassureLipidThickness(grid, sh, d_lipid, true);

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

		vector3 offset(0, 0, 0);

		number h_shift = 2 * (h_corneocyte + d_lipid) / 3;

		// determine max and min x, y dimensions of extruded tkd
		number maxX = 0, maxY = 0, minX = 0, minY = 0;

		for (size_t i = 0; i < shifts.size(); i++) {
			vector3& v = shifts[i].second;
			if (v.x > maxX) {
				maxX = v.x;
			} else if (v.y > maxY) {
				maxY = v.y;
			} else if (v.x < minX) {
				minX = v.x;
			} else if (v.y < minY) {
				minY = v.y;
			}
		}

		number shiftX = maxX + fabs(minX);
		number shiftY = maxY + fabs(minY);

		UG_LOG("min: " << minX << ", " << minY << "\tmax: " << maxX << ", " << maxY << endl);

		for (int i = 0; i < rows; i++) {
			offset.x += shiftX;
			// reset y offset, as we are in a new row
			offset.y = 0;
			// every even row is shifted by 2/3 h_lipid
			if (i % 2) {
				offset.z += h_shift;
			}
			//TODO in every third row, shift z by 1/3 h
			for (int j = 0; j < cols; j++) {
				offset.y += shiftY;
				// reset z offset, as we are in a new col
				offset.z = 0;
				// every even column is shifted by 2/3 h_lipid
				if (j % 2) {
					offset.z += h_shift;
				}

				for (int k = 0; k < high; k++) {
					offset.z += h_corneocyte + d_lipid;
//					UG_LOG(
//							"offset(" << i << ", " << j << ", " << k << "): " << offset << endl);
					Duplicate(grid, sel, offset, aPosition, deselectOld,
							selectNew);
				}
			}
		}

		// delete first created tkd, because it is offset the others due to duplication
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
