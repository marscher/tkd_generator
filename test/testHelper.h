/*
 * testHelper.h
 *
 *  Created on: 03.01.2012
 *      Author: marscher
 */

#ifndef TESTHELPER_H_
#define TESTHELPER_H_

#include "domain_generator.h"
#include "lib_grid/lib_grid.h"
#include <numeric>
using std::vector;

namespace ug {
namespace tkd {

/**
 * calculate volumes of lipid elements
 *
 * @param grid containing tkd with lipid to test
 * @param sh subset handler to assign broken elements to
 * @param d_lipid thickness of lipid matrix to compare with
 * @param assign if true broken elements will be assign to new subset
 */
std::vector<number> meassureLipidVolume(Grid& grid, SubsetHandler& sh,
		number d_lipid, bool assign = false) {
	vector<Face*> faces;
	vector<double> distances;
	Selector sel(grid);

	//	access vertex position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	// select lipid boundary volumes
	SelectAreaBoundary(sel, sh.begin<Volume>(LIPID), sh.end<Volume>(LIPID));

	for (VolumeIterator iter = sh.begin<Volume>(LIPID);
			iter != sh.end<Volume>(LIPID); iter++) {
		Volume* v = *iter;
		distances.push_back(CalculateVolume(v, aaPos));
	}

	return distances;
}

bool checkParallelBoundaryFaces(Grid& grid, SubsetHandler& sh) {
	vector<Face*> faces;
	vector<double> distances;

	//	access vertex position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	Selector sel(grid);
	SelectAreaBoundary(sel, sh.begin<Volume>(LIPID), sh.end<Volume>(LIPID));

	for (VolumeIterator iter = sh.begin<Volume>(LIPID);
			iter != sh.end<Volume>(LIPID); iter++) {
		Volume* v = *iter;
		CollectAssociated(faces, grid, v);
		Face* selFaces[2] = { NULL, NULL };
		for (size_t i_face = 0; i_face < faces.size(); ++i_face) {
			Face* f = faces[i_face];
			if (sel.is_selected(f)) {
				if (selFaces[0]) {
					selFaces[1] = f;
				} else {
					selFaces[0] = f;
				}
			}
		}

		UG_ASSERT(selFaces[0] && selFaces[1],
				"There should be exactly 2 selected faces!");

		vector3 normal1, normal2;
		CalculateNormal(normal1, selFaces[1], aaPos);
		CalculateNormal(normal2, selFaces[0], aaPos);

		VecNormalize(normal1, normal1);
		VecNormalize(normal2, normal2);

		if (!(fabs(VecDot(normal1, normal2) - 1) < SMALL)) {
			return false;
		}
	}
	return true;
}

} // end of namespace tkd
} // end of namespace ug
#endif /* TESTHELPER_H_ */
