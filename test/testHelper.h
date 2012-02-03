/*
 * testHelper.h
 *
 *  Created on: 03.01.2012
 *      Author: marscher
 */

#ifndef TESTHELPER_H_
#define TESTHELPER_H_

#include "../tetrakaidecahedron_generator.h"
#include "../generator.h"
#include "lib_grid/lib_grid.h"

using namespace ug;

namespace tkdGenerator {

number deltaLipidThickness(const vector<number>& distances,
		const number d_lipid) {
	number sum = 0;
	for (size_t i = 0; i < distances.size(); i++) {
		sum += distances[i];
	}
	return fabs(d_lipid / 2 - sum / distances.size());
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
	Selector selDist(grid);

	//	access vertex position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	// select lipid boundary faces
	SelectAreaBoundaryFaces(sel, sh.begin<Volume>(LIPID),
			sh.end<Volume>(LIPID));

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
		UG_LOG("distance : " << dist << endl);
		distances.push_back(dist);
		if (assign && fabs(dist - d_lipid / 2) > 0.1) {
			selDist.select(v);
		}
	}

	return distances;
}

} // end of namespace tkdGenerator
#endif /* TESTHELPER_H_ */
