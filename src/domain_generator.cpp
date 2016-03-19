/*
 * tetrakaidekaeder_generator.cpp
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#include "./domain_generator.h"
#include "util/vecComparator.h"

#include "lib_grid/algorithms/selection_util.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/algorithms/duplicate.h"
#include "lib_grid/grid/grid_object_collection.h"

namespace ug {
namespace tkd {

const number TKDDomainGenerator::REMOVE_DOUBLES_THRESHOLD = 10E-8;

TKDDomainGenerator::TKDDomainGenerator(Domain<3>& d) :
		m_grid(*d.grid()),
		m_sh(*d.subset_handler()),
		m_bSCDomain(true)
{
	setupGridObjects();
}

TKDDomainGenerator::TKDDomainGenerator(Grid& grid, ISubsetHandler& sh) :
		m_grid(grid),
		m_sh(sh),
		m_bSCDomain(true)
{
	if(&grid != sh.grid()) {
		UG_THROW("ERROR: given SubsetHandler not assigned to given Grid instance.");
	}

	setupGridObjects();
}

TKDDomainGenerator::TKDDomainGenerator(Grid& grid, ISubsetHandler& sh,
		bool scDomain) :
		m_grid(grid),
		m_sh(sh),
		m_bSCDomain(scDomain)
{
	if(&grid != sh.grid()) {
		UG_THROW("ERROR: given SubsetHandler not assigned to given Grid instance.");
	}

	setupGridObjects();
}

void TKDDomainGenerator::setGridObject(Grid& grid, ISubsetHandler& sh) {
	if(&grid != sh.grid()) {
		UG_THROW("ERROR: given SubsetHandler not assigned to given Grid instance.");
	}

	this->m_grid = grid;
	this->m_sh = sh;

	setupGridObjects();
}

void TKDDomainGenerator::setupGridObjects() {
	//// set grid options
	m_grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES |
			GRIDOPT_FULL_INTERCONNECTION);
	m_grid.attach_to_vertices(aPosition);
	m_aaPos.access(m_grid, aPosition);
	m_sel.assign_grid(m_grid);
}

void TKDDomainGenerator::setTKDGeometryGenerator(number a, number w, number h,
		bool createLipid, number d_lipid) {
	if (m_pGeomGenerator.get() == NULL)
		m_pGeomGenerator.reset(
				new TKDGeometryGenerator(a, w, h, createLipid, d_lipid));
	else {
		m_pGeomGenerator->setGeometricParams(a, w, h, d_lipid);
	}
}

void TKDDomainGenerator::setSCDomain(bool sc_domain) {
	this->m_bSCDomain = sc_domain;
	if(m_bSCDomain && !getGeometryGenerator().createLipid()) {
		getGeometryGenerator().setCreateLipid(m_bSCDomain);
	}
}

bool TKDDomainGenerator::isSCDomain() const {
	return m_bSCDomain;
}

/**
 * creates lib_grid objects in given grid reference according to given positions and indices
 * @param grid grid instance in which geometric objects will be created
 * @param positions
 * @param indices: numInds1, ind1_1, ind1_2, ..., numInds2, ind2_1, ind2_2, ...
 *
 * the IndexArray contains consecutive integers to describe the geometric objects:
 * First the number of indices of a object is given, followed by the actual indices
 * for the object. Eg. [4, 0, 1, 2, 3] for a tetrahedron.
 * numInds == 4: tetrahedron
 * numInds == 5: pyramid
 * numInds == 6: prism
 * numInds == 8: hexahedron
 */
void TKDDomainGenerator::createGridFromArrays(
		const TKDGeometryGenerator::CoordIndexMap& positions,
		const IndexArray& indices) {
	typedef TKDGeometryGenerator::CoordIndexMap::const_iterator CoordsIter;

	// to lookup created vertices by their index given by geometry generator
	std::map<uint, Vertex*> verts;

	// iterates over all unique vertices and attach position
	for(CoordsIter iter = positions.begin(); iter != positions.end(); ++iter) {
		// create vertex
		Vertex* v = *m_grid.create<RegularVertex>();
		// attach position
		m_aaPos[v] = iter->right;
		// store created vertex by its id
		verts[iter->left] = v;
	}

	// the VolumeDescriptor will be used to create new volumes
	VolumeDescriptor vd;
	// count how much volumes have been created
	uint count = 0;
	// create the elements from the given indices
	m_sh.set_default_subset_index(CORNEOCYTE);
	for (uint i = 0; i < indices.size(); count++) {
		// inner tkd = 63 elements, so lipid matrix volume index begins with 64
		// this assumes that default subset index is set to CORNEOCYTE
		if(count > 62) {
			m_sh.set_default_subset_index(LIPID);
		}

		uint num = indices[i++];
		vd.set_num_vertices(num);

		for (uint j = 0; j < num; ++j) {
			// lookup vertex by its index
			vd.set_vertex(j, verts.at(indices[i++]));
		}

		switch (num) {
		case 4:
			m_grid.create<Tetrahedron>(vd);
			break;
		case 5:
			m_grid.create<Pyramid>(vd);
			break;
		case 6:
			m_grid.create<Prism>(vd);
			break;
		case 8:
			m_grid.create<Hexahedron>(vd);
			break;
		default:
			UG_THROW("wrong number of indices.")
		}
	}

	// assign boundary faces of single tkd
	m_sel.clear();

	// for SC domain assign interface faces of inner tkd to BOUNDARY_CORN and
	// outer faces of lipid volumes to BOUNDARY_LIPID
	if(m_bSCDomain) {
		SelectInterfaceElements(m_sel, m_sh, m_grid.begin<Face>(), m_grid.end<Face>());
		m_sh.assign_subset(m_sel.begin<Face>(), m_sel.end<Face>(), BOUNDARY_CORN);

		GridObjectCollection goc =
				m_sh.get_grid_objects_in_subset(LIPID);
		m_sel.clear();
		SelectBoundaryElements(m_sel, goc.begin<Face>(), goc.end<Face>());
		m_sh.assign_subset(m_sel.begin<Face>(), m_sel.end<Face>(),
				BOUNDARY_LIPID);
	} else {
	// for a simple domain assign outer faces to BOUNDARY_CORN
		GridObjectCollection goc =
						m_sh.get_grid_objects_in_subset(CORNEOCYTE);
		SelectBoundaryElements(m_sel, goc.begin<Face>(), goc.end<Face>());
		m_sh.assign_subset(m_sel.begin<Face>(), m_sel.end<Face>(),
					BOUNDARY_CORN);
	}
}

const TKDDomainGenerator::FaceVec& TKDDomainGenerator::getOppositeFaces(
		const FaceNormalMapping& map,
		const vector3& normal) const {
	vector3 n_flip;
	VecScale(n_flip, normal, -1);
	try {
		return map.at(n_flip);
	} UG_CATCH_THROW("lookup went wrong!")
}

/**
 * assign top/bottom and parallel faces to subset indices.
 * Parallel ones will have consecutive indices and are named after scheme
 * {quad, hex}{a,b}.
 */
void TKDDomainGenerator::assignBoundaryFacesToSubsets(
		const FaceNormalMapping& facesByNormal) {
	const vector3 top(0, 0, 1), bottom(0, 0, -1);

	// start with last used index to assign boundary faces (will be incremented
	// on first matching parallel faces pair)
	int si = BOTTOM;
	// count hexagons and quads for naming
	uint count_hex = 0, count_quad = 0;
	std::stringstream ss;

	for (FNIter iter = facesByNormal.begin(); iter != facesByNormal.end(); ++iter) {
		const vector3& n1 = (*iter).first;
		const FaceVec& faces1 = (*iter).second;
		// calculate normal with flipped orientation and lookup their faces
		const FaceVec& faces2 = getOppositeFaces(facesByNormal, n1);
		UG_ASSERT(faces1.size() == faces2.size(), "face arrays have to match")
		TKDSubsetType t1 = static_cast<TKDSubsetType>(m_sh.get_subset_index(faces1[0])),
			t2 = static_cast<TKDSubsetType>(m_sh.get_subset_index(faces2[0]));

		// ensure faces are not yet assigned (contained in bnd lipid)
		if(t1 >= TOP || t2 >= TOP)
			continue;

		if(n1 == bottom) {
			m_sh.assign_subset(faces1.begin(), faces1.end(), BOTTOM);
			m_sh.assign_subset(faces2.begin(), faces2.end(), TOP);
			AssignAssociatedVerticesToSubset(m_sh, faces1.begin(), faces1.end(), BOTTOM);
			AssignAssociatedVerticesToSubset(m_sh, faces2.begin(), faces2.end(), TOP);

			// todo assign inner edges to subset

			// names of these subsets have to be added here, due to the internal
			// handling of subset handler...
			m_sh.subset_info(TOP).name = "TOP";
			m_sh.subset_info(BOTTOM).name = "BOTTOM";
		} else {
			int a = ++si, b = si + 1;

			if(faces1.size() == 6) {
				ss << "hex" << count_hex++;
			} else {
				ss << "quad" << count_quad++;
			}

			// set subset names
			m_sh.subset_info(a).name = ss.str().append("a");
			m_sh.subset_info(b).name = ss.str().append("b");
			ss.str("");

			// assign faces to subsets
			m_sh.assign_subset(faces1.begin(), faces1.end(), a);
			m_sh.assign_subset(faces2.begin(), faces2.end(), b);
			// and their vertices
			AssignAssociatedVerticesToSubset(m_sh, faces1.begin(), faces1.end(), a);
			AssignAssociatedVerticesToSubset(m_sh, faces2.begin(), faces2.end(), b);

			// todo assign inner edges to subset

			// set subset counter to last subset index used
			si = b;
		}
	}
}

/**
 * calculates the normal for all faces on given subset and stores them in
 * given map.
 * @param shIndex index of outer boundary
 */
void TKDDomainGenerator::mapBoundaryFacesToNormals(
		FaceNormalMapping& facesByNormal, TKDSubsetType shIndex) {
	GridObjectCollection goc = m_sh.get_grid_objects_in_subset(
			shIndex);

	vector3 normal;
	for (FaceIterator iter = goc.begin<Face>(); iter != goc.end<Face>();
			++iter) {
		Face* face = *iter;
		CalculateNormal(normal, face, m_aaPos);
		FaceVec& v = facesByNormal[normal];
		v.push_back(face);
	}
}

/**
 * calculates three vectors for stacking the tkds.
 * Each vector points in the direction defined by the centers of two parallel hexagons.
 */
void TKDDomainGenerator::calculateShiftVectors(UniqueVector3Set& shiftVectors,
		const FaceNormalMapping& facesByNormal) const {
	vector3 c1, c2;

	// find 2 parallel hexagons and calculate vector through their centers,
	// if the vector is unique it will be stored in shiftVectors
	for (FNIter iter = facesByNormal.begin(); iter != facesByNormal.end();
			++iter) {
		const FaceVec& faces = (*iter).second;
		// hexagon?
		if(faces.size() == 6) {
			// get parallel faces
			const FaceVec& parallelHexagon = getOppositeFaces(facesByNormal,
					(*iter).first);
			// calculate their center
			c1 = CalculateCenter(faces.begin(), faces.end(), m_aaPos);
			c2 = CalculateCenter(parallelHexagon.begin(), parallelHexagon.end(),
					m_aaPos);
			// only store unique shifts
			shiftVectors.insert(
					vector3(std::abs(c1.x()) + std::abs(c2.x()),
							std::abs(c1.y()) + std::abs(c2.y()),
							std::abs(c1.z()) + std::abs(c2.z())));
		}
	}
}

/**
 * set subset names
 */
void TKDDomainGenerator::setSubsetNames() {
	m_sh.subset_info(CORNEOCYTE).name = "COR";

	// set subset names only needed if SC domain is wanted
	if(m_bSCDomain) {
		m_sh.subset_info(LIPID).name = "LIP";
		m_sh.subset_info(BOUNDARY_CORN).name = "COR-LIP";
	}
}

/**
 * stacks the created single TKD according to given parameters.
 */
void TKDDomainGenerator::performStacking(uint rows, uint cols, uint layers,
		const FaceNormalMapping& facesByNormal) {
	uint count = rows * cols * layers;
	if(count <= 1)
		return;

	UG_DLOG(MAIN, 0, "creating " << count << " cells." << std::endl);

	m_sel.clear();
	m_sel.enable_autoselection(false);
	// select lipid and corneocyte of first tkd for duplication
	m_sel.select(m_grid.begin<Volume>(), m_grid.end<Volume>());
	m_sel.enable_autoselection(true);

	// do not select duplicates so we can continue to work with initial selection
	bool selectNew = false;
	bool deselectOld = false;

	//// calculate shift vectors for stapeling
	UniqueVector3Set shiftVecs;
	UniqueVector3Set::iterator shiftIter;

	calculateShiftVectors(shiftVecs, facesByNormal);

	vector3 shiftHeight(0, 0, 0), shiftCols(0, 0, 0), shiftRows(0, 0, 0);

	for (shiftIter = shiftVecs.begin(); shiftIter != shiftVecs.end();
			++shiftIter) {
		const vector3& v = *shiftIter;

		if(std::abs(v.x()) < SMALL && std::abs(v.y()) < SMALL) {
			shiftHeight = v;
			shiftHeight.x() = 0;
			shiftHeight.y() = 0;
		} else if(std::abs(v.x()) < SMALL) {
			shiftRows = v;
			shiftRows.x() = 0;
		} else
			shiftCols = v;
	}

	UG_ASSERT(VecLength(shiftHeight) > SMALL && VecLength(shiftCols) > SMALL
			&& VecLength(shiftRows) > SMALL, "shifts not set correctly");

	UG_DLOG(MAIN,0, "shift vectors:\n" << "height: " << shiftHeight
			<< "\tcols: " << shiftCols << "\trows: " << shiftRows << std::endl);

	//// staple in height direction
	vector3 offset_high(0, 0, 0);

	for (uint k = 0; k < layers - 1; k++) {
		VecAdd(offset_high, offset_high, shiftHeight);
		Duplicate(m_grid, m_sel, offset_high, aPosition, deselectOld,
				selectNew);
	}

	//// stack in rows direction
	// select all
	m_sel.select(m_grid.begin<Volume>(), m_grid.end<Volume>());

	vector3 offset_rows(0, 0, 0);

	//	Iterate from 0 to rows - 1 only!
	// every column is shifted by +shiftRows
	// every 3rd column is shifted by -h
	for (uint row = 0; row < rows - 1; row++) {
		VecAdd(offset_rows, offset_rows, shiftRows);
		// next duplicate will complete basis cell
		if((row % 3) == 1) {
			VecSubtract(offset_rows, offset_rows, shiftHeight);
		}
		Duplicate(m_grid, m_sel, offset_rows, aPosition, deselectOld,
				selectNew);
	}

	// remove double vertices for first row
	RemoveDoubles<3>(m_grid, m_grid.vertices_begin(), m_grid.vertices_end(),
			aPosition, REMOVE_DOUBLES_THRESHOLD);

	//// stack in column direction
	// select first generated row
	m_sel.select(m_grid.begin<Volume>(), m_grid.end<Volume>());

	vector3 offset_cols(0, 0, 0);
	vector3 oneThirdHeight;
	VecScale(oneThirdHeight, shiftHeight, 1 / 3.0);
	// loop only to cols - 1, because we already have one column,
	for (uint col = 0; col < cols - 1; col++) {
		VecAdd(offset_cols, offset_cols, shiftCols);

		// shift every second col by -1/3h and -shiftRows.
		if(col % 2) {
			VecSubtract(offset_cols, offset_cols, shiftHeight);
			VecAdd(offset_cols, offset_cols, shiftRows);
		} else {
			VecAdd(offset_cols, offset_cols, oneThirdHeight);
			VecSubtract(offset_cols, offset_cols, shiftRows);
		}
		// in every third row, the y offset has to be reset
		if((col % 3) == 1) {
			offset_cols.y() = 0;
		}

		Duplicate(m_grid, m_sel, offset_cols, aPosition, deselectOld,
				selectNew);
	}

	// finally remove double vertices
	RemoveDoubles<3>(m_grid, m_grid.vertices_begin(), m_grid.vertices_end(),
			aPosition, REMOVE_DOUBLES_THRESHOLD);
}

/**
 * creates one tkd __without__ lipid matrix in given grid reference
 * and stacks it according to given parameters (rows, cols, layers)
 */
void TKDDomainGenerator::createSimpleTKDDomain(number a, number w, number h,
		int rows, int cols, int layers) {
	setSCDomain(false);
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
	UG_DLOG(MAIN, 0, "calling createSCDomain() with following parameter:\n" <<
			"a: " << a << " w: " << w << " h: " << h << " dl: " << d_lipid);
	// check that constraint w > 2a is met
	if((w - 2 * a) <= REMOVE_DOUBLES_THRESHOLD)
		UG_THROW("w > 2a geometric constraint not met!");

	if(std::abs(d_lipid) <= REMOVE_DOUBLES_THRESHOLD)
		UG_THROW("Lipid channel too small.");

	// add call back so dim property of subsets gets updates properly
	m_grid.message_hub()->post_message(
			GridMessage_Creation(GMCT_CREATION_STARTS, 0));

	// deactivate hierarchical insertion for multigrid (creates empty levels...)
	bool hierarchicalInertionEnabled = false;
	MultiGrid* mg = dynamic_cast<MultiGrid*>(&m_grid);
	if(mg) {
		hierarchicalInertionEnabled = mg->hierarchical_insertion_enabled();
		mg->enable_hierarchical_insertion(false);
	}

	setSubsetNames();

	/// create coordinates and vertex indices for 1 tetrakaidecahedron
	if(m_bSCDomain)
		setTKDGeometryGenerator(a, w, h, true, d_lipid);
	else
		setTKDGeometryGenerator(a, w, h, false);

	m_pGeomGenerator->createGeometry();
	/// fill the grid object with coordinates and indices
	createGridFromArrays(m_pGeomGenerator->getPositions(),
						m_pGeomGenerator->getIndices());

	UG_DLOG(MAIN, 0, "Volume of corneocyte: "
			<< m_pGeomGenerator->getVolume(CORNEOCYTE) << std::endl
			<< "volume of lipid: " << m_pGeomGenerator->getVolume(LIPID) << std::endl
			<< "Area of lipid: " << m_pGeomGenerator->getSurface(LIPID) << std::endl);

	TKDSubsetType boundary = BOUNDARY_LIPID;

	if(!m_bSCDomain)
		boundary = BOUNDARY_CORN;


	/// we have to rotate the tkd, because we want to stack the hexagons together
	RotationMatrix r(60);
	TransformVertices(m_grid.begin<Vertex>(), m_grid.end<Vertex>(),
			r, m_aaPos);

	// perform mapping normal -> { faces }
	FaceNormalMapping facesByNormal;
	mapBoundaryFacesToNormals(facesByNormal, boundary);

	/// perform stacking
	performStacking(rows, cols, layers, facesByNormal);

	/// put boundary faces to named subsets according to (quad|hex) and color them
	if (rows * cols * layers == 1)
		assignBoundaryFacesToSubsets(facesByNormal);
	AssignSubsetColors(m_sh);

	// erase temporary subset
//	m_sh.erase_subset(boundary);

	// ensure subsets are proper for simulation
	AdjustSubsetsForSimulation(m_sh, true);

	// activate hierarchical insertion for multigrid again
	if(hierarchicalInertionEnabled && mg)
		mg->enable_hierarchical_insertion(true);

	m_grid.message_hub()->post_message(
					GridMessage_Creation(GMCT_CREATION_STOPS, 0));
}

} // end of namespace tkd
} // end of namespace ug
