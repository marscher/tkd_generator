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

#include <lib_grid/lib_grid.h>

#include "../common_typedefs.h"
#include "../domain_generator.h"
#include "../geometry_generator.h"
#include "../util/volume_calculation.h"

#include <vector>
#include <fstream>

using ug::Grid;
using ug::SubsetHandler;
using ug::Selector;

using std::vector;
using namespace ug::tkd;
// percental allowed deviance of results with calculated values
#define dev_percentage 0.01

struct gridFixture {
	gridFixture() :
					grid(
							VRTOPT_STORE_ASSOCIATED_FACES
									| FACEOPT_STORE_ASSOCIATED_VOLUMES),
					sh(grid) {
		grid.attach_to_vertices(aPosition);
		aaPos = Grid::VertexAttachmentAccessor<APosition>(grid, aPosition);
	}

	Grid grid;
	SubsetHandler sh;
	Grid::VertexAttachmentAccessor<APosition> aaPos;
} domInstance;

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
number calculateSurfaceArea(Grid& grid, SubsetHandler& sh, int subset,
		Grid::VertexAttachmentAccessor<APosition>& aaPos) {

	Selector sel(grid);
	SelectSubsetElements<Face>(sel, sh, subset);

	number sum = 0;

	for (FaceIterator iter = sel.begin<Face>(); iter != sel.end<Face>();
			iter++) {
		Face* f = *iter;

		vector3& a = aaPos[f->vertex(0)];
		vector3& b = aaPos[f->vertex(1)];
		vector3& c = aaPos[f->vertex(2)];
		if (f->num_vertices() == 3) {
			sum += TriangleArea(a, b, c);
		} else if (f->num_vertices() == 4) {
			vector3& d = aaPos[f->vertex(3)];
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
	namespace sp = boost::spirit::classic;
	return sp::parse(str,
	// Begin grammar
			(sp::real_p[sp::push_back_a(v)]
					>> *(',' >> sp::real_p[sp::push_back_a(v)])),
			// End grammar
			sp::space_p).full;
}

void checkVolumeAndSurface(int step, number alpha, number first_vol,
		number first_surface, number firstLipidVol, number a, number h,
		number w) {
	// now check volume and surface calculated with numerically determined results (1 tkd)
	gridFixture gridInstance;
	Grid& grid = gridInstance.grid;
	SubsetHandler& sh = gridInstance.sh;
	Grid::VertexAttachmentAccessor<APosition>& aaPos = gridInstance.aaPos;

	tkd::TKDDomainGenerator gen(grid, sh);
	gen.createSCDomain(a, w, h, 0.1);

	BOOST_REQUIRE(saveToVTK(step,a,h,w));

	// analytical values
	number vol = gen.getGeometryGenerator().getVolume(CORNEOCYTE);
	number lipidVol = gen.getGeometryGenerator().getVolume(LIPID);
	number area = gen.getGeometryGenerator().getSurface(CORNEOCYTE);

	// calculated values
	number area_calc = calculateSurfaceArea(grid, sh, BOUNDARY_CORN, aaPos);

	BOOST_CHECK_CLOSE(area_calc, area, dev_percentage);

	// calc volume over all corneocyte volumes
	number volCalcOnOrignal = CalculateVolume(sh.begin<Volume>(CORNEOCYTE),
			sh.end<Volume>(CORNEOCYTE), aaPos);

	BOOST_CHECK_CLOSE(volCalcOnOrignal, vol, dev_percentage);

	// now tetrahedralize domain and measure all tetrahedron volumes
	eraseVolumes(grid, sh);
	// disable default subset assigning of generated tetrahedrons
	sh.set_default_subset_index(-1);
	std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
	std::cout.rdbuf(NULL);
	BOOST_REQUIRE_MESSAGE(Tetrahedralize(grid, sh, 0, false, false, aPosition),
			"tetrahedralize went wrong");
	// enable outputs again
	std::cout.rdbuf(cout_sbuf);

	// separate lipid and corneocyte subsets again
	SeparateSubsetsByLowerDimSubsets<Volume>(grid, sh);

	int lipidSubset = -1, corneocyteSubset = -1;

	std::vector<Volume*> vols;
	std::vector<Volume*>::iterator voliter;

	// determine subset index of lipid matrix
	for (FaceIterator iter = sh.begin<Face>(BOUNDARY_LIPID);
			iter != sh.end<Face>(BOUNDARY_LIPID); ++iter) {
		Face* f = *iter;
		if (LiesOnBoundary(grid, f)) {
			CollectAssociated(vols, grid, f);
			lipidSubset = sh.get_subset_index(vols[0]);
		}
	}

	if (lipidSubset == 0)
		corneocyteSubset = 1;
	else if (lipidSubset == 1)
		corneocyteSubset = 0;

	number volTet_c = CalculateVolume(sh.begin<Volume>(corneocyteSubset),
			sh.end<Volume>(corneocyteSubset), aaPos);

	number volTet_l = CalculateVolume(sh.begin<Volume>(lipidSubset),
			sh.end<Volume>(lipidSubset), aaPos);

	BOOST_CHECK_CLOSE(volTet_c, vol, dev_percentage);

	// check that tetrahedron volume is close to volume calculated on orignal elements.
	BOOST_CHECK_CLOSE(volTet_c, volCalcOnOrignal, dev_percentage);

	// check that volume / alpha stays constant.
	BOOST_CHECK_CLOSE(volTet_c / alpha, first_vol, dev_percentage);

	// check surface stays constant
	BOOST_CHECK_CLOSE(area_calc, first_surface, dev_percentage);

	// check lipid volume stays constant
	BOOST_CHECK_CLOSE(volTet_l, firstLipidVol, dev_percentage);
	BOOST_CHECK_CLOSE(lipidVol, firstLipidVol, dev_percentage);
}

boost::unit_test::test_suite* initValidateResultsTS(const char* file) {
	boost::unit_test::test_suite* ts = BOOST_TEST_SUITE("validateResults");

	std::ifstream ifs(file);
	if (!ifs.good()) {
		BOOST_MESSAGE("given file: " << file << " not readable.");
		return ts;
	}

	vector<double> v;
	vector<double>::iterator iter;
	char buff[255];

	while (ifs.good()) {
		ifs.getline(buff, 255);
		parse_numbers(buff, v);
	}

	ifs.close();

	uint step = 0;
	number alpha_1_volume = 0, alpha_1_surface = 0, alpha_1_lipid_vol = 0;
	number alpha, a, h, s, w, V, A, V_lip;

	// load swell data from csv and generate geometries
	for (iter = v.begin(); iter != v.end(); step++) {
		UG_LOG("step: " << step << endl);
		alpha = *iter++;
		a = *iter++;
		h = *iter++;
		s = *iter++;
		w = *iter++;
		V = *iter++;
		A = *iter++;
		V_lip = *iter++;
//		UG_LOG("a: " << a << " h: " << h << " s: " <<
//				s << " w: " << "V: " << V << " A: " << A << " Vl: " << V_lip << std::endl);

		assert(a > 0 && s > 0 && w > 0 && h > 0);

		if (step == 0) {
			alpha_1_volume = V;
			alpha_1_surface = A;
			alpha_1_lipid_vol = V_lip;
		}

		boost::unit_test::test_case* tc =
				BOOST_TEST_CASE(boost::bind(checkVolumeAndSurface, step, alpha, alpha_1_volume, alpha_1_surface,
								alpha_1_lipid_vol, a, h, w));

		std::stringstream ss;
		ss << "alpha=" << alpha << " a=" << a << " h=" << h << " w=" << w;
		tc->p_name.set(ss.str());

		ts->add(tc);
	}
	return ts;
}
