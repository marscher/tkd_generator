/*
 * tetrakaidekaeder_generator.h
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#ifndef TETRAKAIDEKAEDER_GENERATOR_H_
#define TETRAKAIDEKAEDER_GENERATOR_H_

#include <vector>

// used ug header
#include "lib_grid/lg_base.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"

#include "registry/registry.h"
#include "common_typedefs.h"
#include "rotation_matrix.h"

namespace tkdGenerator {

using ug::Grid;

/**
 * \param grid
 * \param height
 * \param baseEdgeLength
 * \param diameter
 */
void createGrid(Grid& grid, const CoordsArray& positions,
		const IndexArray& indices);

void TestTKDGenerator(const char* outfile, number height, number baseEdgeLength, number diameter, number d_lipid);

} // end of namespace tkdGenerator

#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */
