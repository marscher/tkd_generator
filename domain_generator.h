/*
 * tetrakaidekaeder_generator.h
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */
#ifndef TETRAKAIDEKAEDER_GENERATOR_H_
#define TETRAKAIDEKAEDER_GENERATOR_H_

namespace tkd {
class TKDDomainGenerator;
class TKDGeometryGenerator;
}

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
	typedef std::set<vector3, vecComperator> shiftSet;
	// access 3d attachment of vertices
	typedef Grid::VertexAttachmentAccessor<APosition> VertexAttachmentAccessor3d;
public:
	TKDDomainGenerator(Grid&, SubsetHandler&);

	void setSubsetHandlerInfo(const char* corneocyte_name,
			const char* lipid_name, const vector4& corneocyte_color,
			const vector4& lipid_color);

	void createTKDDomain(number a, number w, number h, number d_lipid,
			int rows = 1, int cols = 1, int layers = 1);

	TKDGeometryGenerator& getGeometryGenerator() const {
		return *geomGenerator;
	}

	const VertexAttachmentAccessor3d& getVertexAttachmentAccessor() const {
		return aaPos;
	}

private:
	Grid& grid;
	SubsetHandler& sh;
	VertexAttachmentAccessor3d aaPos;
	std::auto_ptr<TKDGeometryGenerator> geomGenerator;
	static const number removeDoublesThreshold = 10E-5;

	void calculateShiftVector(shiftSet& shiftVectors, int subset = LIPID);

	void createGridFromArrays(const CoordsArray& positions,
			const IndexArray& indices);

	void setTKDGeometryGenerator(TKDGeometryGenerator* gen) {
		geomGenerator.reset(gen);
	}
}; // end of class

} // end of namespace tkdGenerator

#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */
