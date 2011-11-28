/*
 * tetrakaidekaeder_generator.h
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#ifndef TETRAKAIDEKAEDER_GENERATOR_H_
#define TETRAKAIDEKAEDER_GENERATOR_H_

#include <vector>

#include <boost/geometry/geometry.hpp>

// used ug header
#include "lib_grid/lg_base.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/file_io/file_io.h"
#include "common/math/ugmath_types.h"
#include "registry/registry.h"

namespace trans = boost::geometry::strategy::transform;
using boost::geometry::dsv;

namespace tkdGenerator {

typedef boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian> point_type;

typedef std::vector<ug::vector3> CoordsArray;
typedef std::vector<int> IndexArray;
typedef const ug::vector3& vec3Ref;
typedef double number;

/**
 * \param grid
 * \param height
 * \param baseEdgeLength
 * \param diameter
 */
void GenerateTetrakaidecahedron(ug::Grid& grid, number& height,
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

void createPrism(vec3Ref v1, vec3Ref v2, vec3Ref v3,
				 vec3Ref v4, vec3Ref v5, vec3Ref v6,
				 CoordsArray& posOut, IndexArray& indsOut);

void createTetrahedron(vec3Ref v1, vec3Ref v2, vec3Ref v3, vec3Ref v4,
		CoordsArray& posOut, IndexArray& indsOut);

void TestTKDGenerator(const char* outfile, number height, number baseEdgeLength, number diameter);

} // end of namespace tkdGenerator

#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */
