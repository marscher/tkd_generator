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

#include <map>

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
	uint count = 1;
	Volume* vol = NULL;

	// create the elements from the given indices
	for (size_t i = 0; i < indices.size(); count++) {
		int num = indices[i++];
		vd.set_num_vertices(num);
		for (int j = 0; j < num; ++j)
			vd.set_vertex(j, vertices[indices[i++]]);

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

		// inner tkd = 63 elements, so lipid matrix volume index begins with 64
		// this assumes that default subset index is set to CORNEOCYTE
		if (count >= 64) {
			sh.assign_subset(vol, LIPID);
			for (uint k = 0; k < vol->num_faces(); k++) {
				Face* f = grid.get_face(vol, k);
				sh.assign_subset(f, LIPID);
			}
		}
	}

	// remove double vertices
	RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
			aPosition, 10E-5);
}

/**
 * calculates three vectors for stacking the tkds.
 * Each vector points in the direction defined by the centers of two parallel hexagons.
 */
void calculateShiftVector(shiftSet& shiftVectors, Grid& grid, SubsetHandler& sh,
		Grid::VertexAttachmentAccessor<APosition>& aaPos, int subset) {
	// collect faces associated to unique normal
	map<vector3, vector<Face*>, vecComperator> facesByNormal;
	map<vector3, vector<Face*>, vecComperator>::iterator fnIter;
	vector3 normal, c1, c2;

	for (FaceIterator iter = sh.begin<Face>(subset);
			iter != sh.end<Face>(subset); iter++) {
		Face* face = *iter;
		if (IsBoundaryFace3D(grid, face)) {
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
			vector<Face*>& parallelHexagon = facesByNormal[inverseNormal];
			// calculate their center
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

/**
 * set Subset informations
 * @param sh SubsetHandler reference to fill with given information
 * @param corneocyte_name default "corneocytes"
 * @param lipid_name default "lipid matrix"
 * @param corneocyte_color default blue
 * @param lipid_color default green
 */
void setSubsetHandlerInfo(SubsetHandler& sh, const char* corneocyte_name,
		const char* lipid_name, const vector4& corneocyte_color,
		const vector4& lipid_color) {

	// Subset 1 for corneocytes (same as in Feuchters tkdmodeller)
	sh.set_default_subset_index(CORNEOCYTE);
	// fixme crashs subsethandler
	sh.subset_info(CORNEOCYTE).name = corneocyte_name;
	sh.subset_info(LIPID).name = lipid_name;

	sh.subset_info(LIPID).color = lipid_color;
	sh.subset_info(CORNEOCYTE).color = corneocyte_color;
}

/**
 * creates one tkd with lipid matrix in given grid reference
 * and stacks it according to given parameters (rows, cols, layers)
 */
void createTKDDomain(Grid& grid, SubsetHandler& sh, number baseEdgeLength,
		number diameter, number height, number d_lipid, int rows,
		int cols, int layers) {

	UG_LOG(
			"a: " << baseEdgeLength << " w: " << diameter << " h: " << height << " d: " << d_lipid << endl);

	//// set grid options
	grid.set_options(GRIDOPT_STANDARD_INTERCONNECTION);
	grid.attach_to_vertices(aPosition);

	//// create coordinates and vertex indices for 1 tetrakaidecahedron
	TKDGeometryGenerator generator(height, baseEdgeLength, diameter,
			d_lipid);
	generator.createDomain();
	//// fill the grid object
	createGridFromArrays(grid, sh, generator.getPositions(),
			generator.getIndices());

	//// perform stacking
	uint count = rows * cols * layers;
	if (count > 1) {
		UG_LOG("creating " << count << " cells with lipid matrix." << endl);

		Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

		Selector sel(grid);
		sel.enable_autoselection(false);
		// select lipid and corneocyte of first tkd for duplication
		sel.select(grid.begin<Volume>(), grid.end<Volume>());

		// do not select duplicates so we can continue to work with initial selection
		bool selectNew = false;
		bool deselectOld = false;

		//// calculate shift vectors for stapeling
		shiftSet shiftVecs;
		shiftSet::iterator shiftIter;

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

		UG_LOG(
				"shiftHeight: " << shiftHeight << "\tshiftcols: " << shiftCols << "\tshiftrows: " << shiftRows << endl);

		UG_ASSERT(
				VecLength(shiftHeight)> SMALL && VecLength(shiftCols)> SMALL && VecLength(shiftRows)> SMALL,
				"shifts not set correctly");

		//// staple in height direction
		vector3 offset_high(0, 0, 0);

		for (int k = 0; k < layers - 1; k++) {
			VecAdd(offset_high, offset_high, shiftHeight);
			Duplicate(grid, sel, offset_high, aPosition, deselectOld,
					selectNew);
		}

		//// stack in rows direction
		// select all
		sel.select(grid.begin<Volume>(), grid.end<Volume>());

		vector3 offset_rows(0, 0, 0);

		//	Iterate from 0 to rows - 1 only!
		// every column is shifted by +shiftRows
		// every 3rd column is shifted by -h
		for (int row = 0; row < rows - 1; row++) {
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

		//// stack in column direction
		// select first generated row
		sel.select(grid.begin<Volume>(), grid.end<Volume>());

		vector3 offset_cols(0, 0, 0);
		vector3 oneThirdHeight;
		VecScale(oneThirdHeight, shiftHeight, 1 / 3.0);
		// loop only to cols - 1, because we already have one column,
		for (int col = 0; col < cols - 1; col++) {
			VecAdd(offset_cols, offset_cols, shiftCols);

			// shift every second col by -1/3h and -shiftRows.
			if (col % 2) {
				VecSubtract(offset_cols, offset_cols, shiftHeight);
				VecAdd(offset_cols, offset_cols, shiftRows);
			} else {
				VecAdd(offset_cols, offset_cols, oneThirdHeight);
				VecSubtract(offset_cols, offset_cols, shiftRows);
			}
			// in every third row, the y offset has to be reset
			if ((col % 3) == 1) {
				offset_cols.y = 0;
			}

			Duplicate(grid, sel, offset_cols, aPosition, deselectOld,
					selectNew);
		}

		// finally remove double vertices
		RemoveDoubles<3>(grid, grid.vertices_begin(), grid.vertices_end(),
				aPosition, 10E-5);
	}
}

/**
 * fills given grid reference with tkd domain with given parameters.
 * Provides a SubsetHandler with default parameters.
 */
SubsetHandler createTKDDomainDefaultSubsetInfo(Grid& grid, number height,
		number baseEdgeLength, number diameter, number d_lipid, int rows,
		int cols, int layers) {

	SubsetHandler sh(grid);

	const char* corneocyte_name = "corneocytes";
	const char* lipid_name = "lipid matrix";
	const vector4 corneocyte_color = vector4(0, 1, 0, 0);
	const vector4 lipid_color = vector4(0, 0, 1, 0);

	setSubsetHandlerInfo(sh, corneocyte_name, lipid_name, corneocyte_color,
			lipid_color);

	createTKDDomain(grid, sh, baseEdgeLength, diameter, height, d_lipid, rows,
			cols, layers);

	return sh;
}

// register tkd generator functions for usage in ug_script
extern "C" void InitUGPlugin(ug::bridge::Registry* reg,
		std::string parentGroup) {
	std::string grp(parentGroup);
	grp.append("tkd_generator/");

	//	add TKD-Generator method
	reg->add_function(
			"CreateTKDDomain",
			&createTKDDomainDefaultSubsetInfo,
			grp /*,
			"returns a SubsetHandler and a Grid with stacked Tetrakaidecahedrons." */);

	reg->add_function("CreateTKDDomain_", &createTKDDomain, grp/*, "Grid with stacked Tetrakaidecahedrons."*/);

	reg->add_function("InitSubsetHandler", &setSubsetHandlerInfo, grp);
}

} // end of namespace tkdGenerator
