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

// fixme params should be const, but can't because they're calling non const workers
number CalculateVolume(const Tetrahedron&);
number CalculateVolume(const Prism&);
number CalculateVolume(const Pyramid&);
number CalculateVolume(const Hexahedron&);
number CalculateVolume(const Volume&);
number CalculateVolume(VolumeIterator& begin, VolumeIterator& end);
} // end namespace ug

// include implementation
//#include "volume_calculation_impl.h"

#endif /* VOLUME_CALCULATION_H_ */
