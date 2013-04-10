/*
 * tetrakaidekaeder_generator.h
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */
#ifndef TETRAKAIDEKAEDER_GENERATOR_H_
#define TETRAKAIDEKAEDER_GENERATOR_H_

#include "lib_grid/grid/grid.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/tools/selector_grid.h"
#include "lib_disc/domain.h"

#include "common_typedefs.h"
#include "geometry_generator.h"
// for std::less<ug::MathVector<T, dim> >
#include "util/vecComparator.h"

#include <map>
#include <set>

namespace ug {
namespace tkd {

/**
 *
 */
class TKDDomainGenerator {

public:
	// unique vector set
	typedef std::set<vector3> UniqueVector3Set;
	// access 3d attachment of vertices
	typedef Grid::VertexAttachmentAccessor<APosition> VertexAttachmentAccessor3d;

	// store faces associated to their normal
	typedef std::vector<Face*> FaceVec;
	// this uses std::less<vector3> to sort keys (see util/vecComparator.h)
	typedef std::map<vector3, FaceVec> FaceNormalMapping;
	typedef FaceNormalMapping::const_iterator FNIter;

	TKDDomainGenerator(Domain<3>& domain3d);

	/**
	 *
	 */
	TKDDomainGenerator(Grid&, ISubsetHandler&);

	/**
	 *
	 */
	TKDDomainGenerator(Grid&, ISubsetHandler&, bool scDomain,
			bool distinctBndSubsetInds = true);

	/**
	 * sets grid and subset handler to use
	 */
	void setGridObject(Grid&, ISubsetHandler&);

	/**
	 * sets whether a stratum corneum (SC) domain should be created (nested tkd's).
	 */
	void setSCDomain(bool);

	/**
	 * will this domain generator create a SC domain?
	 */
	bool isSCDomain() const;

	void createSCDomain(number a, number w, number h, number d_lipid, int rows =
			1, int cols = 1, int layers = 1);

	void createSimpleTKDDomain(number a, number w, number h, int rows = 1,
			int cols = 1, int layers = 1);

	VertexAttachmentAccessor3d& getVertexAttachmentAccessor() {
		return m_aaPos;
	}

	TKDGeometryGenerator& getGeometryGenerator() {
		if(m_pGeomGenerator.get() == NULL)
			m_pGeomGenerator.reset(new TKDGeometryGenerator());
		return *m_pGeomGenerator;
	}

	/// helper functions for registry to get subset indices
	static int getLipidIndex() {
		return LIPID;
	}

	static int getCorneocyteIndex() {
		return CORNEOCYTE;
	}

	static int getCorneocyteBoundaryIndex() {
		return BOUNDARY_CORN;
	}

	static int getLipidBoundaryIndex() {
		return BOUNDARY_LIPID;
	}

private:
	Grid& m_grid;
	ISubsetHandler& m_sh;
	Selector m_sel;
	VertexAttachmentAccessor3d m_aaPos;
	std::auto_ptr<TKDGeometryGenerator> m_pGeomGenerator;
	/**
	 * threshold for lib_grid's <RemoveDoubles>() to indicate when objects
	 * are too close to each other and will be removed.
	 */
	static const number REMOVE_DOUBLES_THRESHOLD;
	/**
	 * indicates whether a stratum corneum domain is beeing created (set in ctor)
	 */
	bool b_scDomain;

	/**
	 * returns the faces container in given map with opposite of given normal
	 */
	inline const FaceVec& getOppositeFaces(const FaceNormalMapping& map,
			const vector3& normal) const;

	/**
	 * associate faces of given tkd subset type to their normal
	 */
	void mapBoundaryFacesToNormals(FaceNormalMapping&,
			TKDSubsetType bnd = BOUNDARY_LIPID);

	/**
	 * will assign two consequtive (n, n+1) subset indices for each pair
	 * of parallel boundary faces to use them e.g. for periodic identification
	 */
	void assignBoundaryFacesToSubsets(const FaceNormalMapping&);

	void calculateShiftVectors(UniqueVector3Set&, const FaceNormalMapping&) const;

	void createGridFromArrays(const TKDGeometryGenerator::CoordIndexMap&,
			const IndexArray&);

	/**
	 * Sets the geometric parameters of the encapsulated TKDGeometryGenerator
	 * to given params.
	 * By default a lipid matrix with a thickness of 0.1 will be created.
	 * @param a base edge lenth of base hexagon of tkd
	 * @param w 'diameter' of tkd
	 * @param h height of tkd
	 * @param createLipid indicate whether a lipid matrix should be created.
	 * @param d_lipid thickness of lipid matrix
	 */
	void setTKDGeometryGenerator(number a, number w, number h,
			bool createLipid = true, number d_lipid = 0.1);
	/**
	 * set up common grid objects
	 */
	void setupGridObjects();

	void setSubsetNames();
};

} // end of namespace tkd
} // end of namespace ug
#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */
