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

using ug::vector3;
using ug::Grid;

/**
 * \param grid
 * \param height
 * \param baseEdgeLength
 * \param diameter
 */
void GenerateTetrakaidecahedron(Grid& grid, number& height,
		number& baseEdgeLength, number& diameter);

/**
 * indsOut: numInds1, ind1_1, ind1_2, ..., numInds2, ind2_1, ind2_2, ...
 *
 * numInds == 4: tetrahedron
 * numInds == 5: pyramid
 * numInds == 6: prism
 * numInds == 8: hexahedron
 */
void GenerateTetrakaidecahedron(CoordsArray&, IndexArray&,
							  number& height, number& baseEdgeLength, number& diameter);


void TestTKDGenerator(const char* outfile, number height, number baseEdgeLength, number diameter);

} // end of namespace tkdGenerator

#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */
