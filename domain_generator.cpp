/*
 * tetrakaidekaeder_generator.cpp
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#include "domain_generator.h"
#include "test/testHelper.h"

#include "lib_grid/lib_grid.h"

#include <map>
namespace ug {
namespace tkd {

const number TKDDomainGenerator::removeDoublesThreshold = 10E-5;

TKDDomainGenerator::TKDDomainGenerator(Grid& grid, ISubsetHandler& sh) :
		m_grid(grid), m_sh(sh), b_scDomain(true) {
	if (&grid != sh.grid()) {
		UG_THROW(
				"ERROR: given SubsetHandler not assigned to given Grid instance.");
	}

	m_grid.attach_to_vertices(aPosition);
	m_aaPos.access(m_grid, aPosition);

	m_grid.set_options(GRIDOPT_AUTOGENERATE_SIDES);
}

TKDDomainGenerator::TKDDomainGenerator(Grid& grid, ISubsetHandler& sh, bool scDomain) :
		m_grid(grid), m_sh(sh), b_scDomain(scDomain) {
	if (&grid != sh.grid()) {
		UG_THROW(
				"ERROR: given SubsetHandler not assigned to given Grid instance.");
	}

	m_grid.attach_to_vertices(aPosition);
	m_aaPos.access(m_grid, aPosition);

	m_grid.set_options(GRIDOPT_AUTOGENERATE_SIDES);
}

void TKDDomainGenerator::setGridObject(Grid& grid, ISubsetHandler& sh) {
	if (&grid != sh.grid()) {
		UG_THROW(
				"ERROR: given SubsetHandler not assigned to given Grid instance.");
	}

	this->m_grid = grid;
	this->m_sh = sh;
	this->m_grid.attach_to_vertices(aPosition);
	m_aaPos.access(m_grid, aPosition);
	m_grid.set_options(GRIDOPT_AUTOGENERATE_SIDES);
}

/**
 * sets whether a stratum corneum domain should be created (nested tkd's).
 */
void TKDDomainGenerator::setIsSCDomain(bool sc_domain) {
	this->b_scDomain = sc_domain;
	if(b_scDomain && !getGeometryGenerator().createLipid()) {
		getGeometryGenerator().setCreateLipid(b_scDomain);
	}
}

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
void TKDDomainGenerator::createGridFromArrays(const CoordsArray& positions,
		const IndexArray& indices, bool sc_domain) {
	// generate vertices in the grid and store them in an array, so that we can index them
	std::vector<VertexBase*> vertices(positions.size());
	for (size_t i = 0; i < positions.size(); ++i) {
		VertexBase* v = *m_grid.create<Vertex>();
		m_aaPos[v] = positions[i];
		vertices[i] = v;
	}
	// the VolumeDescriptor will be used to create new volumes
	VolumeDescriptor vd;
	uint count = 0;
	// create the elements from the given indices
	for (size_t i = 0; i < indices.size(); count++) {
		// inner tkd = 63 elements, so lipid matrix volume index begins with 64
		// this assumes that default subset index is set to CORNEOCYTE
		if (count > 62) {
			m_sh.set_default_subset_index(LIPID);
		}
		int num = indices[i++];
		vd.set_num_vertices(num);
		for (int j = 0; j < num; ++j)
			vd.set_vertex(j, vertices[indices[i++]]);
		switch (num) {
		case 4:
			m_grid.create<ug::Tetrahedron>(vd);
			break;
		case 5:
			m_grid.create<ug::Pyramid>(vd);
			break;
		case 6:
			m_grid.create<ug::Prism>(vd);
			break;
		case 8:
			m_grid.create<ug::Hexahedron>(vd);
			break;
		}
	}

	// remove double vertices
	RemoveDoubles<3>(m_grid, m_grid.vertices_begin(), m_grid.vertices_end(),
			aPosition, removeDoublesThreshold);

	m_sh.assign_subset(m_grid.begin<VertexBase>(), m_grid.end<VertexBase>(), -1);
	m_sh.assign_subset(m_grid.begin<EdgeBase>(), m_grid.end<EdgeBase>(), -1);
	m_sh.assign_subset(m_grid.begin<Face>(), m_grid.end<Face>(), -1);

	// assign boundary faces
	if(sc_domain) {
		Selector sel(m_grid);
		SelectInterfaceElements(sel, m_sh, m_grid.begin<Face>(), m_grid.end<Face>());
		m_sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), BOUNDARY_CORN);

		sel.clear();
		GeometricObjectCollection goc = m_sh.get_geometric_objects_in_subset(LIPID);
		SelectBoundaryElements(sel, goc.begin<Face>(), goc.end<Face>());
		m_sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), BOUNDARY_LIPID);
	}

	// TODO: assign TOP/BOTTOM
	/*if(false) {
		// TOP/BOT
		Selector sel(m_grid);
		vector3 normal;

		GeometricObjectCollection goc = m_sh.get_geometric_objects_in_subset(BOUNDARY_LIPID);
		for (FaceIterator iter = goc.begin<Face>();
			  iter != goc.end<Face>(); iter++) {
					Face* face = *iter;
					CalculateNormal(normal, face, m_aaPos);
					double norm=sqrt(normal*normal);
					double dir=normal*vector3(0,0,1.0);
					dir /=norm;
					if (dir==norm) {}
			}
	}*/
	AdjustSubsetsForSimulation(m_sh, true);
}

/**
 * calculates three vectors for stacking the tkds.
 * Each vector points in the direction defined by the centers of two parallel hexagons.
 */
void TKDDomainGenerator::calculateShiftVector(shiftSet& shiftVectors,
		TKDSubsetType shIndex) {
	// collect faces associated to unique normal
	map<vector3, vector<Face*>, vecComparator> facesByNormal;
	map<vector3, vector<Face*>, vecComparator>::iterator fnIter;
	vector3 normal, c1, c2;

	GeometricObjectCollection goc = m_sh.get_geometric_objects_in_subset(shIndex);
	for (FaceIterator iter = goc.begin<Face>();
			iter != goc.end<Face>(); iter++) {
		Face* face = *iter;
			CalculateNormal(normal, face, m_aaPos);
			facesByNormal[normal].push_back(face);
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
			c1 = CalculateCenter(faces.begin(), faces.end(), m_aaPos);
			c2 = CalculateCenter(parallelHexagon.begin(), parallelHexagon.end(),
					m_aaPos);
			// only store unique shifts
			using std::abs;
			shiftVectors.insert(vector3(abs(c1.x) + abs(c2.x),
					abs(c1.y) + abs(c2.y),
					abs(c1.z) + abs(c2.z)));
		}
	}
}

/**
 * set Subset informations
 * @param corneocyte_name
 * @param lipid_name default
 * @param corneocyte_color
 * @param lipid_color
 */
void TKDDomainGenerator::setSubsetHandlerInfo(const char* corneocyte_name,
		const char* lipid_name, const vector4& corneocyte_color,
		const vector4& lipid_color) {

	// Subset 1 for corneocytes (same as in Feuchters tkdmodeller)
	m_sh.set_default_subset_index(CORNEOCYTE);
	m_sh.subset_info(CORNEOCYTE).name = corneocyte_name;
	m_sh.subset_info(CORNEOCYTE).color = corneocyte_color;

	m_sh.subset_info(LIPID).name = lipid_name;
	m_sh.subset_info(LIPID).color = lipid_color;

	// set subset boundary infos
	vector4 boundary_color;
	VecMultiply(boundary_color, lipid_color, 0.6);
	m_sh.subset_info(BOUNDARY_CORN).name = "COR-LIP";
	m_sh.subset_info(BOUNDARY_CORN).color = boundary_color;

	VecMultiply(boundary_color, corneocyte_color, 0.6);
	m_sh.subset_info(BOUNDARY_LIPID).name = "TOP";
	m_sh.subset_info(BOUNDARY_LIPID).color = boundary_color;

}

/**
 * creates one tkd __without__ lipid matrix in given grid reference
 * and stacks it according to given parameters (rows, cols, layers)
 */
void TKDDomainGenerator::createSimpleTKDDomain(number a, number w, number h,
		int rows, int cols, int layers) {
	setIsSCDomain(false);
	createSCDomain(a, w, h, -1, rows, cols, layers);
}

/**
 * creates one tkd __with__ lipid matrix in given grid reference
 * and stacks it according to given parameters (rows, cols, layers)
 * @param a
 * @param w
 * @param h
 * @param d_lipid
 * @param rows
 * @param cols
 * @parm layers
 */
void TKDDomainGenerator::createSCDomain(number a, number w, number h,
		number d_lipid, int rows, int cols, int layers) {




	UG_LOG("calling createSCDomain() with following parameter:\n" <<
			"a: " << a << " w: " << w << " h: " << h << " dl: " << d_lipid << endl);
	// check that constraint w > 2a is met
	if (! std::abs(w - 2 * a) > removeDoublesThreshold)
		UG_THROW("w > 2a geometric constraint not met!");


	// deactivate hierarchical insertion for multigrid (creates empty levels...)
	bool hierarchicalInertionEnabled = false;
	MultiGrid* mg= dynamic_cast<MultiGrid*>(&m_grid);
	if (mg)
	{
		hierarchicalInertionEnabled = mg->hierarchical_insertion_enabled();
		mg->enable_hierarchical_insertion(false);
	}


	if(b_scDomain)
		setSubsetHandlerInfo("corneocytes", "lipid matrix", vector4(0, 1, 0, 0),
				vector4(0, 0, 1, 0));

	//// set grid options
	m_grid.set_options(GRIDOPT_STANDARD_INTERCONNECTION);
	m_grid.attach_to_vertices(aPosition);
	//// create coordinates and vertex indices for 1 tetrakaidecahedron
	if(b_scDomain)
		setTKDGeometryGenerator(a, w, h, true, d_lipid);
	else
		setTKDGeometryGenerator(a, w, h, false);

	geomGenerator->createGeometry();
	//// fill the grid object with coordinates and indices
	createGridFromArrays(geomGenerator->getPositions(),
			geomGenerator->getIndices(), b_scDomain);
	UG_DLOG(APP, 1,
			"Volume of corneocyte: " << geomGenerator->getVolume(CORNEOCYTE) << endl
			<< "volume of lipid: " << geomGenerator->getVolume(LIPID) << endl
			<< "Area of lipid: " << geomGenerator->getSurface(LIPID) << endl);
	//// perform stacking
	uint count = rows * cols * layers;

	if (count > 1) {
		UG_DLOG(APP, 1,
				"creating " << count << " cells with lipid matrix." << endl);

		Grid::VertexAttachmentAccessor<APosition> m_aaPos(m_grid, aPosition);

		Selector sel(m_grid);
		sel.enable_autoselection(false);
		// select lipid and corneocyte of first tkd for duplication
		sel.select(m_grid.begin<Volume>(), m_grid.end<Volume>());

		// do not select duplicates so we can continue to work with initial selection
		bool selectNew = false;
		bool deselectOld = false;

		//// calculate shift vectors for stapeling
		shiftSet shiftVecs;
		shiftSet::iterator shiftIter;
		TKDSubsetType boundary = BOUNDARY_LIPID;

		if(!b_scDomain)
			boundary = BOUNDARY_CORN;

		calculateShiftVector(shiftVecs, boundary);

		vector3 shiftHeight(0, 0, 0), shiftCols(0, 0, 0), shiftRows(0, 0, 0);

		for (shiftIter = shiftVecs.begin(); shiftIter != shiftVecs.end();
				shiftIter++) {
			using std::abs;
			const vector3& v = *shiftIter;

			if (abs(v.x) < SMALL && abs(v.y) < SMALL) {
				shiftHeight = v;
				shiftHeight.x = 0;
				shiftHeight.y = 0;
			} else if (abs(v.x) < SMALL) {
				shiftRows = v;
				shiftRows.x = 0;
			} else
				shiftCols = v;
		}

		UG_ASSERT(
				VecLength(shiftHeight) > SMALL && VecLength(shiftCols) > SMALL && VecLength(shiftRows) > SMALL,
				"shifts not set correctly");

		UG_DLOG(APP, 1,
				"shift vectors:\n" << "height: " << shiftHeight << "\tcols: "
				<< shiftCols << "\trows: " << shiftRows << endl)

		//// staple in height direction
		vector3 offset_high(0, 0, 0);

		for (int k = 0; k < layers - 1; k++) {
			VecAdd(offset_high, offset_high, shiftHeight);
			Duplicate(m_grid, sel, offset_high, aPosition, deselectOld,
					selectNew);
		}

		//// stack in rows direction
		// select all
		sel.select(m_grid.begin<Volume>(), m_grid.end<Volume>());

		vector3 offset_rows(0, 0, 0);

		//	Iterate from 0 to rows - 1 only!
		// every column is shifted by +shiftRows
		// every 3rd column is shifted by -h
		for (int row = 0; row < rows - 1; row++) {
			VecAdd(offset_rows, offset_rows, shiftRows);
			if ((row % 3) == 1) {
				VecSubtract(offset_rows, offset_rows, shiftHeight);
			}

			Duplicate(m_grid, sel, offset_rows, aPosition, deselectOld,
					selectNew);
		}

		// remove double vertices for first row
		RemoveDoubles<3>(m_grid, m_grid.vertices_begin(), m_grid.vertices_end(),
				aPosition, removeDoublesThreshold);

		//// stack in column direction
		// select first generated row
		sel.select(m_grid.begin<Volume>(), m_grid.end<Volume>());

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

			Duplicate(m_grid, sel, offset_cols, aPosition, deselectOld,
					selectNew);
		}

		// finally remove double vertices
		RemoveDoubles<3>(m_grid, m_grid.vertices_begin(), m_grid.vertices_end(),
				aPosition, removeDoublesThreshold);
	}



	// activate hierarchical insertion for multigrid again
	if(hierarchicalInertionEnabled && mg)
		mg->enable_hierarchical_insertion(true);

	// TODO: grid _changed_calback??
	UpdateMaxDimensionOfSubset(m_sh, "dim");
}


} // end of namespace tkdGenerator
} // end of namespace ug
