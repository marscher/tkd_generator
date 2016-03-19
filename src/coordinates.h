#ifndef __coordinates_h__
#define __coordinates_h__

#include "common_typedefs.h"

namespace ug {
namespace tkd {

/// \addtogroup tkd_generator
/// \{

void CalculateLipidCoords(CoordsArray& l, double a_corneo, double high,
		double width, double d_lipid, const vector3& offset);

// end group tkd_generator
/// \}

} // end of namespace tkd
} // end of namespace ug
#endif
