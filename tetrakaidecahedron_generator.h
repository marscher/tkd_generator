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

namespace tkdGenerator {

/**
 * \param grid
 * \param sh
 * \param height
 * \param baseEdgeLength
 * \param diameter
 */
void createGridFromArrays(ug::Grid& grid, ug::SubsetHandler& sh, const CoordsArray& positions,
		const IndexArray& indices);

void TestTKDGenerator(const char* outfile, number height, number baseEdgeLength,
		number diameter, number d_lipid, int rows, int cols, int high);

} // end of namespace tkdGenerator

#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */
