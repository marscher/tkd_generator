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
#include <set>

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

		//// calculate shift vectors for stapeling
		set<vector3, vecComperator> shiftVecs;
		set<vector3, vecComperator>::iterator shiftIter;

		calculateShiftVector(shiftVecs, grid, sh, aaPos);
		vector3 shiftHeight(0, 0, 0), shiftCols(0, 0, 0), shiftRows(0, 0, 0);

		for (shiftIter = shiftVecs.begin(); shiftIter != shiftVecs.end();
				shiftIter++) {
			const vector3& v = *shiftIter;

			if (fabs(v.x) < SMALL && fabs(v.y) < SMALL) {
				shiftHeight = v;
				shiftHeight.x = 0;
				shiftHeight.y = 0;
			} else if (fabs(v.x) < SMALL) {
				shiftRows = v;
				shiftRows.x = 0;
			} else
				shiftCols = v;
		}
/*
// Das hier ist nicht gut. Warum erstellst du eine Kante?
// Das fŸhrt natŸrlich zum Absturz, weil du eine Kante ohne vertices erstellst.
// Dann weist du denen aber was zu... Auf Mac crashed das schon bei der zuweisung.
// Lass das doch einfach raus.
		EdgeBase* e = *grid.create<Edge>();
		aaPos[e->vertex(1)] = origin;
		aaPos[e->vertex(1)] = shiftCols;
*/
// this crashes on ug::GridWriterUGX::add_elements_to_node(rapidxml::xml_node<char>*, ug::Grid&) only if subset is assigned
//		sh.assign_subset(e, 3);

		UG_ASSERT(
				VecLength(shiftHeight)> SMALL && VecLength(shiftCols)> SMALL && VecLength(shiftRows)> SMALL,
				"shifts not set correctly");

		UG_LOG(
				"shiftHeight: " << shiftHeight << "\tshiftcols: " << shiftCols << "\tshiftrows: " << shiftRows << endl);

		//// staple in height direction
		vector3 offset_high(0, 0, 0);

		for (uint k = 0; k < high - 1; k++) {
			VecAdd(offset_high, offset_high, shiftHeight);
			Duplicate(grid, sel, offset_high, aPosition, deselectOld,
					selectNew);
		}

//Delete is not necessary. Iterate from 0 to high - 1 instead!
		// delete first created tkd, because it is offset the others due to duplication
		//grid.erase(sel.begin<Vertex>(), sel.end<Vertex>());
		//grid.erase(sel.begin<Edge>(), sel.end<Edge>());
		//grid.erase(sel.begin<Volume>(), sel.end<Volume>());

		//// staple in rows direction
		// select all
		sel.select(grid.begin<Volume>(), grid.end<Volume>());

		vector3 offset_rows(0, 0, 0);

//	Iterate from 0 to rows - 1 only!
		// every column is shifted by +shiftRows
		// every 3rd column is shifted by -h
		for (uint row = 0; row < rows - 1; row++) {
			VecAdd(offset_rows, offset_rows, shiftRows);
			if ((row % 3) == 1) {
				VecSubtract(offset_rows, offset_rows, shiftHeight);
			}

			Duplicate(grid, sel, offset_rows, aPosition, deselectOld,
					selectNew);
		}

		// remove double vertices for first row
		RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
				aPosition, 10E-5);

		//// staple in column direction
		// select first generated row
		sel.select(grid.begin<Volume>(), grid.end<Volume>());

		vector3 offset_cols(0, 0, 0);
		vector3 oneThirdHeight, twoThirdHeight;
		VecScale(oneThirdHeight, shiftHeight, 1 / 3.0);
		VecScale(twoThirdHeight, shiftHeight, 2 / 3.0);
		// loop only to cols - 1, because we already have one column,
		// so first col has to be shifted by -shiftRows, +1/3 h and +shiftCols
		// second by -1/3 h and ...
		for (uint col = 0; col < cols - 1; col++) {
			UG_LOG("col: " << col << endl);
			VecAdd(offset_cols, offset_cols, shiftCols);
			VecAdd(offset_cols, offset_cols, oneThirdHeight);
			VecSubtract(offset_cols, offset_cols, shiftRows);

			// fixme shift every third col by -1/3h and -shiftrows.
			if (col % 2) {
				UG_LOG("sub 1/3h " << endl);
				VecSubtract(offset_cols, offset_cols, oneThirdHeight);
//				VecSubtract(offset_cols, offset_cols, shiftRows);
			}
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
