/*
 * validateSwellResults.cpp
 *
 *  Created on: Mar 8, 2012
 *      Author: marscher
 */

#include <boost/test/unit_test.hpp>
// parser
#include <boost/spirit/include/classic.hpp>

#include "../domain_generator.h"
#include "../geometry_generator.h"
#include "../util/volume_calculation.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

//using namespace boost;
using namespace boost::spirit::classic;
using namespace std;
using namespace ug;

/**
 * delete all volumes and keep surface areas of the 2 nested tkds,
 * so we can perform a tetrahedralisation.
 */
void eraseVolumes(Grid& grid, SubsetHandler&sh) {
	Selector sel(grid);
	SelectBoundaryElements(sel, grid.begin<Face>(), grid.end<Face>());
	SelectAssociatedGeometricObjects(sel);

	for (VolumeIterator viter = sh.begin<Volume>(LIPID);
			viter != sh.end<Volume>(LIPID); viter++) {
		Volume* vol = *viter;
		for (size_t i_face = 0; i_face < vol->num_faces(); i_face++) {
			Face* f = grid.get_face(vol, i_face);
			bool selectedVertices = false;
			for (size_t i = 0; i < f->num_vertices(); i++) {
				if (sel.is_selected(f->vertex(i))) {
					selectedVertices = true;
					break;
				}
			}
			if (!selectedVertices) {
				sel.select(f);
			}
		}
	}

	SelectAssociatedGeometricObjects(sel);
	InvertSelection(sel);
	EraseSelectedObjects(sel);
}

bool parse_numbers(char const* str, vector<double>& v) {
	return parse(str,
	// Begin grammar
			(real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)])),
			// End grammar
			space_p).full;
}

BOOST_AUTO_TEST_CASE(validateSwellResults) {
	vector<double> v;
	vector<double>::iterator iter;

	ifstream ifs("/home/marscher/workspace/ug4/data/sols_a10_w30_h1_.csv");
	char buff[255];
	while (ifs.good()) {
		ifs.getline(buff, 255);
		parse_numbers(buff, v);
	}
	ifs.close();
	number dev_percentage = 3;

	// load swell data from csv and generate geometries
	for (iter = v.begin(); iter != v.end();) {
		double alpha = *iter++;
		double a = *iter++;
		double h = *iter++;
		double w = *iter++;
		double d = 0.1;

		if (a < 0 || w < 0 || w < 2 * a || h < 0)
			continue;

//		alpha: 10 a: 1 h: 30 w: 10.0025a: 1 w: 10.0025 h: 30 dl: 0.1
		//fixme führt zu a = 0 im tkd generator ...

		UG_LOG(
				"alpha: " << alpha << " a: " << a << " h: " << h << " w: " << w);

		Grid grid;
		SubsetHandler sh(grid);
		tkd::TKDDomainGenerator gen(grid, sh);
		gen.createTKDDomain(a,w,h,d);

		number vol = gen.getGeometryGenerator().getVolume(CORNEOCYTE);

		Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

		// calc volume over all corneocyte volumes
		number volCalcOnOrignal = CalculateVolume(sh.begin<Volume>(CORNEOCYTE),
				sh.end<Volume>(CORNEOCYTE), aaPos);
		// allow 1% deviance
		BOOST_CHECK_CLOSE(volCalcOnOrignal, vol, dev_percentage);

		// now tetrahedralize domain and meassure all tetrahedron volumes
		eraseVolumes(grid, sh);
		// disable default subset assigning of generated tetrahedrons
		sh.set_default_subset_index(-1);
		BOOST_REQUIRE_MESSAGE(
				Tetrahedralize(grid, sh, 0, false, false, aPosition),
				"tetrahedralize went wrong");
		// separate lipid and corneocyte subsets again
		SeparateSubsetsByLowerDimSubsets<Volume>(grid, sh);

		// determine which subset index is lipid
		int lipidSubset = -1;
		for (FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>();
				++iter) {
			Face* f = *iter;
			if (LiesOnBoundary(grid, f)) {
				lipidSubset = sh.get_subset_index(f);
				break;
			}
		}

//		UG_LOG(
//				"lip index: " << lipidSubset << " name: "<< sh.get_subset_name(lipidSubset) << endl);
//		UG_LOG("name_c: " << sh.get_subset_name(lipidSubset+1) << endl);

// subset index of corneocyte = lipid subset + 1
		number volTet = CalculateVolume(sh.begin<Volume>(lipidSubset + 1),
				sh.end<Volume>(lipidSubset + 1), aaPos);
		BOOST_CHECK_CLOSE(volTet, vol, dev_percentage);

		// check that tetrahedron volume is close to volume calculated on orignal elements.
		BOOST_CHECK_CLOSE(volTet, volCalcOnOrignal, dev_percentage);

//		stringstream s;
//		s << "/tmp/tkd/tkd_" << alpha << ".ugx";
//		SaveGridToFile(grid, sh, s.str().c_str());
//		break;
	}
}
