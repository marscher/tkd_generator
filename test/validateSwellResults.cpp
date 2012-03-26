/*
 * validateSwellResults.cpp
 *
 *  Created on: Mar 8, 2012
 *      Author: marscher
 */
#include "validateSwellResults.h"

#include <boost/test/unit_test.hpp>
// parser
#include <boost/spirit/include/classic.hpp>

#include "../util/VTKOutputGrid.h"

#include "../domain_generator.h"
#include "../geometry_generator.h"
#include "../util/volume_calculation.h"

#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace boost::spirit::classic;
using namespace std;
using namespace ug;

struct gridFixture {
	gridFixture() :
			grid(VRTOPT_STORE_ASSOCIATED_FACES), sh(grid), gen(grid, sh) {
		grid.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(grid, aPosition);
	}

	Grid grid;
	SubsetHandler sh;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
	tkd::TKDDomainGenerator gen;
};

number quadArea(const vector3& a, const vector3& b, const vector3& c,
		const vector3&d) {
	vector3 ab, ac, cross;
	VecSubtract(ab, b, a);
	VecSubtract(ac, c, a);
	VecCross(cross, ab, ac);
	return VecLength(cross);
}

/**
 * calculates surface area of lipid = sum of triangles + sum of quads
 * uses function quadArea.
 */
number calculateSurfaceArea(Grid& grid, SubsetHandler& sh,
		Grid::VertexAttachmentAccessor<APosition>& aaPos) {
	BOOST_CHECKPOINT("begin calc surface");
	Selector sel(grid);

	SelectBoundaryElements(sel, grid.begin<Face>(), grid.end<Face>());

	number sum = 0;

	for (FaceIterator iter = sel.begin<Face>(); iter != sel.end<Face>();
			iter++) {
		Face* f = *iter;

		vector3 a, b, c;
		a = aaPos[f->vertex(0)];
		b = aaPos[f->vertex(1)];
		c = aaPos[f->vertex(2)];
		if (f->num_vertices() == 3) {
			sum += TriangleArea(a, b, c);
		} else if (f->num_vertices() == 4) {
			vector3 d;
			d = aaPos[f->vertex(3)];
			sum += quadArea(a, b, c, d);
		}
	}
	return sum;
}

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

void checkResults(int step, number a, number h, number w) {
	gridFixture fix;
	Grid& grid = fix.grid;
	SubsetHandler& sh = fix.sh;
	Grid::VertexAttachmentAccessor<APosition>& aaPos = fix.aaPos;
	tkd::TKDDomainGenerator& gen = fix.gen;
	gen.createTKDDomain(a, w, h, 0.1);

	int dev_percentage = 1;

	stringstream stream;
//	stream << "/tmp/tkd/tkd_" << alpha << ".vtk";
//	SaveGridToFile(grid, sh, stream.str().c_str());
//	print_subset(&grid, sh, "/tmp/test.vtk", 0,1,0,true);
	BOOST_REQUIRE(SaveGridToVTK(grid, sh, "/tmp/tkd/testvtk", aaPos, step));

	number vol = gen.getGeometryGenerator().getVolume(CORNEOCYTE);
	number area = gen.getGeometryGenerator().getSurface(LIPID);

	number area_calc = calculateSurfaceArea(grid, sh, aaPos);

	BOOST_CHECK_CLOSE(area_calc, area, dev_percentage);

// calc volume over all corneocyte volumes
	number volCalcOnOrignal = CalculateVolume(sh.begin<Volume>(CORNEOCYTE),
			sh.end<Volume>(CORNEOCYTE), aaPos);
//		UG_LOG("volume of corneocyte (org): " << volCalcOnOrignal << endl)

	BOOST_CHECK_CLOSE(volCalcOnOrignal, vol, dev_percentage);

// now tetrahedralize domain and meassure all tetrahedron volumes
	eraseVolumes(grid, sh);
// disable default subset assigning of generated tetrahedrons
	sh.set_default_subset_index(-1);
	BOOST_REQUIRE_MESSAGE( Tetrahedralize(grid, sh, 0, false, false, aPosition),
			"tetrahedralize went wrong");
// separate lipid and corneocyte subsets again
	SeparateSubsetsByLowerDimSubsets<Volume>(grid, sh);

	int lipidSubset = -1, corneocyteSubset = -1;
	for (FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>();
			++iter) {
		Face* f = *iter;
		if (LiesOnBoundary(grid, f)) {
			lipidSubset = sh.get_subset_index(f);
//				UG_LOG("found boundary face" << endl);
			sh.assign_subset(f, 2);
			break;
		}
	}

	if (lipidSubset == 0)
		corneocyteSubset = lipidSubset + 1;
	else if (lipidSubset == 1)
		corneocyteSubset = lipidSubset - 1;

	number volTet_c = CalculateVolume(sh.begin<Volume>(corneocyteSubset),
			sh.end<Volume>(corneocyteSubset), aaPos);

	number volTet_l = CalculateVolume(sh.begin<Volume>(lipidSubset),
			sh.end<Volume>(lipidSubset), aaPos);

	if (volTet_c < volTet_l) {
		UG_LOG(
				"swapping lipid and corneo volume of tetrahedralisation!" << endl)
		swap(volTet_c, volTet_l);
	}

//		UG_LOG("volume of corneocyte(tet): " << volTet_c << endl)
	BOOST_CHECK_CLOSE(volTet_c, vol, 1);

// check that tetrahedron volume is close to volume calculated on orignal elements.
	BOOST_CHECK_CLOSE(volTet_c, volCalcOnOrignal, dev_percentage);
}

boost::unit_test::test_suite* validateResultsTS(const char* file) {
	boost::unit_test::test_suite* ts = BOOST_TEST_SUITE("validateResults");

	ifstream ifs(file);
	if (!ifs.good()) {
		BOOST_MESSAGE("file not readable.");
	}

	vector<double> v;
	vector<double>::iterator iter;
	char buff[255];

	while (ifs.good()) {
		ifs.getline(buff, 255);
		parse_numbers(buff, v);
	}
	ifs.close();
	int step = 0;
// load swell data from csv and generate geometries
	for (iter = v.begin(); iter != v.end();step++) {
		double alpha = *iter++;
		double a = *iter++;
		double h = *iter++;
		double t = *iter++;
		double asl = *iter++;

		number s = h / 3 * t;
		number w = (2 * sqrt(3) * a + 3 * s) / sqrt(3);
		if (a < 0 || s < 0 || w < 0)
			continue;
		ts->add(BOOST_TEST_CASE(boost::bind(checkResults, step, a, h, w)));
	}
	return ts;
}
