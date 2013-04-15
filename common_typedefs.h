/*
 * common_typedefs.h
 *
 *  Created on: 06.12.2011
 *      Author: marscher
 */

#ifndef COMMON_TYPEDEFS_H_
#define COMMON_TYPEDEFS_H_

#include <vector>
#include <common/types.h>
#include "common/math/ugmath_types.h"
#include "common/math/misc/math_util.h"

namespace ug {
namespace tkd {

/// \addtogroup tkd_generator
/// \{

typedef std::vector<vector3> CoordsArray;
typedef std::vector<int> IndexArray;

/// subset indices in created domains subset handler
enum TKDSubsetType {
	/// index for lipid volumes (volumes surrounding inner tkd)
	LIPID = 0,
	/// index for corneocyte (inner tkd)
	CORNEOCYTE,
	/// boundary faces of corneocyte (inner tkd)
	BOUNDARY_CORN,
	/// boundary faces of lipid volumes
	BOUNDARY_LIPID,
	/// top hexagon (in z direction)
	TOP,
	/// bottom hexagon (in z direction)
	BOTTOM
};

// end group tkd_generator
/// \}

} // end of namespace tkd
} // end of namespace ug
#endif /* COMMON_TYPEDEFS_H_ */
