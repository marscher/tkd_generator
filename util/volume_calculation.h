/*
 * volume_calculation.h
 *
 *  Created on: 11.01.2012
 *      Author: marscher
 */

#ifndef VOLUME_CALCULATION_H_
#define VOLUME_CALCULATION_H_

#include "lib_grid/lg_base.h"

namespace ug {

number CalculateVolume(const Tetrahedron&);
number CalculateVolume(const Prism&);
number CalculateVolume(const Pyramid&);
number CalculateVolume(const Hexahedron&);
number CalculateVolume(const Volume&);
number CalculateVolume(geometry_traits<Volume>::iterator, geometry_traits<Volume>::iterator);

} // end namespace ug

#endif /* VOLUME_CALCULATION_H_ */
