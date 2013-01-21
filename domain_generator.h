/*
 * tetrakaidekaeder_generator.h
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */
#ifndef TETRAKAIDEKAEDER_GENERATOR_H_
#define TETRAKAIDEKAEDER_GENERATOR_H_

#include "lib_grid/lib_grid.h"
#include "common_typedefs.h"
#include "geometry_generator.h"

#include <set>

namespace ug {
namespace tkd {

class TKDDomainGenerator {

	// unique vector set
	typedef std::set<vector3> UniqueVector3Set;
	// access 3d attachment of vertices
	typedef Grid::VertexAttachmentAccessor<APosition> VertexAttachmentAccessor3d;

public:
	TKDDomainGenerator(Grid&, ISubsetHandler&);
	TKDDomainGenerator(Grid&, ISubsetHandler&, bool scDomain);

	void setGridObject(Grid&, ISubsetHandler&);

	void setIsSCDomain(bool);

	void setSubsetHandlerInfo(const char* corneocyte_name,
			const char* lipid_name, const vector4& corneocyte_color,
			const vector4& lipid_color);

	void createSCDomain(number a, number w, number h, number d_lipid, int rows =
			1, int cols = 1, int layers = 1);

	void createSimpleTKDDomain(number a, number w, number h, int rows = 1,
			int cols = 1, int layers = 1);

	const VertexAttachmentAccessor3d& getVertexAttachmentAccessor() const {
		return m_aaPos;
	}

	TKDGeometryGenerator& getGeometryGenerator() {
		if(geomGenerator.get() == NULL)
			geomGenerator.reset(new TKDGeometryGenerator());
		return *geomGenerator;
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
	std::auto_ptr<TKDGeometryGenerator> geomGenerator;
	static const number REMOVE_DOUBLES_THRESHOLD;
	bool b_scDomain;

	void assignPeriodicBoundaryFacesToSubsets(std::map<vector3, std::vector<Face*> >&);

	void calculateShiftVectors(UniqueVector3Set& shiftVectors,
			std::map<vector3, std::vector<Face*> >& facesByNormal, TKDSubsetType sh =
			BOUNDARY_LIPID);

	void createGridFromArrays(const CoordsArray& positions,
			const IndexArray& indices, bool sc_domain);

	void setTKDGeometryGenerator(number a, number w, number h,
			bool createLipid = true, number d_lipid = -1) {
		if(geomGenerator.get() == NULL)
			geomGenerator.reset(
					new TKDGeometryGenerator(a, w, h, createLipid, d_lipid));
		else {
			geomGenerator->setGeometricParams(a, w, h, d_lipid);
		}
	}
	void setupGridObjects();
};

} // end of namespace tkd
} // end of namespace ug
#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */
