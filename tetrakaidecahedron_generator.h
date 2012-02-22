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
#include <set>

namespace tkdGenerator {

using namespace ug;
using namespace std;

// unique vector set
typedef std::set<vector3, vecComperator> shiftSet;

// subset indices
enum Subsets {
	LIPID = 0, CORNEOCYTE
};

void createGridFromArrays(Grid& grid, SubsetHandler& sh,
		const CoordsArray& positions, const IndexArray& indices);

void calculateShiftVector(shiftSet& shiftVectors, Grid& grid, SubsetHandler& sh,
		Grid::VertexAttachmentAccessor<APosition>& aaPos, int subset = LIPID);

void setSubsetHandlerInfo(SubsetHandler& sh, const char* corneocyte_name,
		const char* lipid_name, const vector4& corneocyte_color,
		const vector4& lipid_color);

void createTKDDomain(Grid& grid, SubsetHandler& sh, number baseEdgeLength,
		number diameter, number height, number d_lipid, int rows, int cols,
		int layers);

SubsetHandler createTKDDomainDefaultSubsetInfo(Grid& grid, number height,
		number baseEdgeLength, number diameter, number d_lipid, int rows,
		int cols, int layers);

} // end of namespace tkdGenerator

#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */
