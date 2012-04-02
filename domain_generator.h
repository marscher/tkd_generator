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
#include "util/vecComparator.h"
#include "geometry_generator.h"

#include <set>

namespace tkd {

using namespace ug;
using namespace std;

class TKDDomainGenerator {

	// unique vector set
	typedef std::set<vector3, vecComparator> shiftSet;
	// access 3d attachment of vertices
	typedef Grid::VertexAttachmentAccessor<APosition> VertexAttachmentAccessor3d;

public:
	TKDDomainGenerator(Grid&, SubsetHandler&);

	void setGridObject(Grid&, SubsetHandler&);

	void setSubsetHandlerInfo(const char* corneocyte_name,
			const char* lipid_name, const vector4& corneocyte_color,
			const vector4& lipid_color);

	void createTKDDomain(number a, number w, number h, number d_lipid,
			int rows = 1, int cols = 1, int layers = 1);

	VertexAttachmentAccessor3d getVertexAttachmentAccessor() const {
		return aaPos;
	}

	tkd::TKDGeometryGenerator& getGeometryGenerator() {
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
	Grid& grid;
	SubsetHandler& sh;
	VertexAttachmentAccessor3d aaPos;
	std::auto_ptr<TKDGeometryGenerator> geomGenerator;
	static const number removeDoublesThreshold;

	void calculateShiftVector(shiftSet& shiftVectors);

	void createGridFromArrays(const CoordsArray& positions,
			const IndexArray& indices);

	void setTKDGeometryGenerator(number a, number w, number h, number d_lipid) {
		if (geomGenerator.get() == NULL)
			geomGenerator.reset(new TKDGeometryGenerator(a, w, h, d_lipid));
		else {
			geomGenerator->setGeometricParams(a,w,h,d_lipid);
		}
	}
}; // end of class

}// end of namespace tkdGenerator

#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */
