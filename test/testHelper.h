/*
 * testHelper.h
 *
 *  Created on: 03.01.2012
 *      Author: marscher
 */

#ifndef TESTHELPER_H_
#define TESTHELPER_H_

#include "../domain_generator.h"
#include "lib_grid/lib_grid.h"
#include <numeric>

using namespace ug;

namespace ug {
namespace tkd {

number deltaLipidThickness(const vector<number>& distances,
		const number d_lipid) {
	number sum = std::accumulate(distances.begin(), distances.end(), 0);
	if (sum > 10E-6)
		return fabs(d_lipid / 2 - sum / distances.size());
	return 0;
}

/**
 * @param grid containing tkd with lipid to test
 * @param sh subset handler to assign broken elements to
 * @param d_lipid thickness of lipid matrix to compare with
 * @param assign if true broken elements will be assign to new subset
 */
vector<number> meassureLipidThickness(Grid& grid, SubsetHandler& sh,
		number d_lipid, bool assign = false) {
	vector<Face*> faces;
	vector<number> distances;
	Selector sel(grid);

	//	access vertex position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	// select lipid boundary volumes
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

		vector3 center1 = CalculateCenter(selFaces[0], aaPos);
		vector3 center2 = CalculateCenter(selFaces[1], aaPos);
		vector3 normal;
		CalculateNormal(normal, selFaces[1], aaPos);

		vector3 point;
		ProjectPointToPlane(point, center1, center2, normal);
		number dist = VecDistance(center1, point);
		distances.push_back(dist);
		// mark pairs of wrong distanced faces
		if (assign && fabs(dist - d_lipid / 2) > 0.1) {
			sh.assign_subset(selFaces[0], 3);
			sh.assign_subset(selFaces[1], 3);
		}
	}

	return distances;
}

bool checkParallelBoundaryFaces(Grid& grid, SubsetHandler& sh) {
	vector<Face*> faces;
	vector<number> distances;

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
