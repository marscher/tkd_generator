/*
 * tetrakaidekaeder_generator.cpp
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#include "tetrakaidecahedron_generator.h"
#include "generator.h"
#include "test/testHelper.h"
#include "util/vecComparator.h"
#include "lib_grid/lib_grid.h"
#include "registry/registry.h"

#include <map>

namespace tkdGenerator {

using namespace ug;
using namespace std;

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

	sh.assign_subset(grid.begin<Vertex>(), grid.end<Vertex>(), CORNEOCYTE);
	sh.assign_subset(grid.begin<Face>(), grid.end<Face>(), CORNEOCYTE);
	sh.assign_subset(grid.begin<Volume>(), grid.end<Volume>(), CORNEOCYTE);
}

void generateLipidMatrixForSingleTKD(Grid& grid, SubsetHandler& sh,
		number h_scale, Grid::VertexAttachmentAccessor<APosition>& aaPos) {

	vector<Face*> faces;

	// select boundaryFaces for extrusion
	for (FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>();
			++iter) {
		if (IsVolumeBoundaryFace(grid, *iter))
			faces.push_back(*iter);
	}

	// select vertices which are beeing extruded
	Selector sel(grid);
	sel.enable_autoselection(true);

	// extruding boundary faces
	Extrude(grid, NULL, NULL, &faces, origin);

	// assign all extruded geometric objects to LIPID subset
	sh.assign_subset(sel.begin<Vertex>(), sel.end<Vertex>(), LIPID);
	sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), LIPID);
	sh.assign_subset(sel.begin<Volume>(), sel.end<Volume>(), LIPID);

	// scale every extruded vertices with factor h_scale
	for (VertexBaseIterator iter = sel.begin<VertexBase>();
			iter != sel.end<VertexBase>(); iter++) {
		VertexBase* v = *iter;
		vector3& pos = aaPos[v];
		VecScale(pos, pos, h_scale);
	}
}

// should work for both extruded and non extruded single tkd centered around z axis.
void calculateShiftVector(set<vector3, vecComperator>& shiftVectors, Grid& grid,
SubsetHandler& sh, Grid::VertexAttachmentAccessor<APosition>& aaPos,
int subset = LIPID) {
	// collect faces associated to unique normal
	map<vector3, std::vector<Face*>, vecComperator> facesByNormal;
	map<vector3, vector<Face*>, vecComperator>::iterator fnIter;

	for (FaceIterator iter = sh.begin<Face>(subset);
			iter != sh.end<Face>(subset); iter++) {
		Face* face = *iter;

		if (IsBoundaryFace3D(grid, face)) {
			vector3 normal;
			CalculateNormal(normal, face, aaPos);
			facesByNormal[normal].push_back(face);
		}
	}

	// find 2 parallel hexagons and calculate vector through their centers,
	// if the vector is unique it will be stored in shiftVectors
	for (fnIter = facesByNormal.begin(); fnIter != facesByNormal.end();
			fnIter++) {
		vector<Face*>& faces = (*fnIter).second;

		// hexagon?
		if (faces.size() == 6) {
			vector3 normal = (*fnIter).first;
			vector3 inverseNormal = normal;
			// switch orientation
			VecScale(inverseNormal, inverseNormal, -1);
			// get parallel faces
			vector<Face*> parallelHexagon = facesByNormal[inverseNormal];
			// calculate their center
			vector3 c1, c2;
			c1 = CalculateCenter(faces.begin(), faces.end(), aaPos);
			c2 = CalculateCenter(parallelHexagon.begin(), parallelHexagon.end(),
					aaPos);
			// only store unique shifts
			vector3 tmp(fabs(c1.x) + fabs(c2.x), fabs(c1.y) + fabs(c2.y),
					fabs(c1.z) + fabs(c2.z));
			shiftVectors.insert(tmp);
		}
	}
}

void logfkt(const number& d) {
	UG_LOG("d: " << d << endl)
}

void logfkt_vec(const vector3& v) {
	UG_LOG("v: " << v << endl)
}

void GenerateCorneocyteWithLipid(Grid& grid, SubsetHandler& sh,
		number a_corneocyte, number w_corneocyte, number h_corneocyte,
		number d_lipid, uint rows, uint cols, uint high) {

	UG_LOG(
			"a: " << a_corneocyte << " w: " << w_corneocyte << " h: " << h_corneocyte << " d: " << d_lipid << endl);

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

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	// rotate tkd by 60 degrees so side hexagons are orientated parallel to y axis
	RotationMatrix R(60);
	TransformVertices(grid.begin<Vertex>(), grid.end<Vertex>(), R, aaPos);

	number h_scale = (h_corneocyte + d_lipid) / h_corneocyte;

	// generate lipid matrix
	generateLipidMatrixForSingleTKD(grid, sh, h_scale, aaPos);

//	vector<number> dist = meassureLipidThickness(grid, sh, d_lipid, false);
//	for_each(dist.begin(), dist.end(), logfkt);

	UG_ASSERT(checkParallelBoundaryFaces(grid, sh),
			"corresponding faces are not parallel!");

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

		number h = (h_corneocyte + d_lipid);
		number h23 = 2. / 3. * h;
		number h3 = 1. / 3. * h;

		set<vector3, vecComperator> shiftVecs;
		set<vector3, vecComperator>::iterator shiftIter;

		calculateShiftVector(shiftVecs, grid, sh, aaPos);
		vector3 shiftHeight(0, 0, 0), shiftCols(0, 0, 0), shiftRows(0, 0, 0);

		for (shiftIter = shiftVecs.begin(); shiftIter != shiftVecs.end();
				shiftIter++) {
			const vector3& v = *shiftIter;

			if (fabs(v.x) < SMALL && fabs(v.y < SMALL))
				shiftHeight = v;
			else if (fabs(v.x) < SMALL)
				shiftRows = v;
			else
				shiftCols = v;
		}

		UG_ASSERT(
				VecLength(shiftHeight)> SMALL && VecLength(shiftCols)> SMALL && VecLength(shiftRows)> SMALL,
				"shifts not set correctly");

		for_each(shiftVecs.begin(), shiftVecs.end(), logfkt_vec);

		UG_LOG(
				"shiftHeight: " << shiftHeight << "\tshiftcols: " << shiftCols << "\tshiftrows: " << shiftRows << endl);

		uint countRow = 1;
		number z = 0;
		vector3 offset_rows(0, 0, 0);

		for (uint row = 0; row < rows; row++, countRow++) {
			z = -h;

			switch (countRow) {
			// every second row is shifted by 1/3 h
			case 2:
				z += h3;
				break;
				// very third row is shifted by 2/3 h
			case 3:
				z += h23;
				countRow = 0;
				break;
			}

			VecScale(offset_rows, shiftRows, row + 1);

			offset_rows.z = z;

			for (uint k = 0; k < high; k++) {
				offset_rows.z += h;
				Duplicate(grid, sel, offset_rows, aPosition, deselectOld,
						selectNew);
			}
		}

		// delete first created tkd, because it is offset the others due to duplication
		grid.erase(sel.begin<Vertex>(), sel.end<Vertex>());
		grid.erase(sel.begin<Edge>(), sel.end<Edge>());
		grid.erase(sel.begin<Volume>(), sel.end<Volume>());

		// remove double vertices for first row
		RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
				aPosition, 10E-5);

		// select first generated row
		sel.select(grid.begin<Volume>(), grid.end<Volume>());

		vector3 offset_cols(0, 0, 0);
		for (uint col = 0; col < cols - 1; col++) {

			// every second col is shifted by 1/3 h and -1/2 shiftY
			if ((col % 2) == 0) {
				z = h3; // correct
			} else {
				z = 0;
			}

			VecScale(offset_cols, shiftCols, col+1);
			offset_cols.z += z;

			Duplicate(grid, sel, offset_cols, aPosition, deselectOld,
					selectNew);
		}

		// finally remove double vertices
		RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
				aPosition, 10E-5);
	}
}

////////////////////////////////////////////////////////////////////////
///	test tetrakaidekahedron generator
// fixme upstream bug in registry, recursive template instantiation if rows etc. is given as uint!
void TestTKDGenerator(const char *outfile, number height, number baseEdgeLength,
		number diameter, number d_lipid, int rows, int cols, int high) {

	Grid g;
	SubsetHandler sh;

	GenerateCorneocyteWithLipid(g, sh, baseEdgeLength, diameter, height,
			d_lipid, rows, cols, high);

	SaveGridToFile(g, sh, outfile);

	UG_LOG("vec comparisions: " << vecComperator::comps << endl)
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
