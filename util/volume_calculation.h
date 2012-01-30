/*
 * volume_calculation.h
 *
 *  Created on: 11.01.2012
 *      Author: marscher
 */

#ifndef VOLUME_CALCULATION_H_
#define VOLUME_CALCULATION_H_

#include "lib_grid/lib_grid.h"

namespace ug {

number CalculateVolume(const Tetrahedron&, Grid::VertexAttachmentAccessor<APosition>&);
number CalculateVolume(const Prism&, Grid::VertexAttachmentAccessor<APosition>&);
number CalculateVolume(const Pyramid&, Grid::VertexAttachmentAccessor<APosition>&);
number CalculateVolume(const Hexahedron&, Grid::VertexAttachmentAccessor<APosition>&);
number CalculateVolume(const Volume&, Grid::VertexAttachmentAccessor<APosition>&);
number CalculateVolume(geometry_traits<Volume>::iterator,
					   geometry_traits<Volume>::iterator,
					   Grid::VertexAttachmentAccessor<APosition>&);

} // end namespace ug

#endif /* VOLUME_CALCULATION_H_ */
