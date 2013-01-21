/*
 * testLipidDistance.cpp
 *
 *  Created on: 10.01.2012
 *      Author: marscher
 */
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "common/util/path_provider.h"
#include "common/util/file_util.h"

#include "testHelper.h"
#include "../domain_generator.h"

#include <utility>
#include <cmath>

using namespace boost::unit_test;
//using namespace tkd;
using namespace ug;
using namespace ug::tkd;

BOOST_AUTO_TEST_SUITE(lipid)

/**
 * performs checks if lipid thickness (= distance all boundary faces to inner faces) is correct
 */
void checkThickness(number a, number h, number w, number d_lipid) {
	// init libgrid objects
	Grid grid(GRIDOPT_STANDARD_INTERCONNECTION);
	grid.attach_to_vertices(aPosition);
	SubsetHandler sh(grid);

	// create grid with new tkd generator
	TKDDomainGenerator gen(grid, sh);
	gen.createSCDomain(a,w,h,d_lipid);

	// calculate distances
	vector<number> dist = meassureLipidThickness(grid, sh, d_lipid);
	number delta = deltaLipidThickness(dist, d_lipid);

//	BOOST_MESSAGE(
//			"max dist: " << *std::max_element(dist.begin(), dist.end())
//			<< " min dist: " << *std::min_element(dist.begin(), dist.end()) << "\n");
//
//	BOOST_MESSAGE(
//			"abweichung dlipid new: " << delta << "\n");
//
//	// check if average abbreviation of lipid thickness is below given threshold
//	BOOST_MESSAGE("testing if deltaLipidThickness(a=" << a << ", h="
//			<< h << ", w=" << w << ", d=" << d_lipid << ") is small.");
	BOOST_CHECK_SMALL(delta, 10E-6);
}

BOOST_AUTO_TEST_CASE(compareThicknessWithOld) {
	number a = 10;
	number h = sqrt(6) * a;
	number w = 3 * a;
	number d_lipid = 8;

	// init libgrid objects
	Grid gridNew(GRIDOPT_STANDARD_INTERCONNECTION);
	Grid gridOld(GRIDOPT_STANDARD_INTERCONNECTION);
	SubsetHandler shold(gridOld);
	SubsetHandler shnew(gridNew);

	// create grid with new tkd generator
	gridNew.attach_to_vertices(aPosition);
	gridOld.attach_to_vertices(aPosition);

	TKDDomainGenerator gen(gridNew, shnew);
	gen.createSCDomain(a,w,h,d_lipid);

	// grid created with same params with old tkdmodeller
	string oldfile =
	PathProvider::get_path(ROOT_PATH)
	+ "/../plugins/experimental/tkd_generator/test/nucleus3d_10-24.495-30-8.000_1x1x1.ugx";

	bool gridLoaded = LoadGridFromFile(gridOld, shold, oldfile.c_str());
	BOOST_REQUIRE_MESSAGE(gridLoaded,
			"old grid file " << oldfile << " could not be loaded.");

	BOOST_MESSAGE("old tkdmodeler: ________________________________");
	vector<number> distold = meassureLipidThickness(gridOld, shold, d_lipid);
	BOOST_MESSAGE("tkdgeometrygenerator: ________________________________");
	vector<number> distnew = meassureLipidThickness(gridNew, shnew, d_lipid);
	BOOST_MESSAGE("old tkdmodeler: ________________________________");
	BOOST_MESSAGE(
			"max dist: " << *std::max_element(distold.begin(), distold.end()) << " min dist: " << *std::min_element(distold.begin(), distold.end()) << "\n");
	BOOST_MESSAGE("tkdgeometrygenerator: ________________________________");
	BOOST_MESSAGE(
			"max dist: " << *std::max_element(distnew.begin(), distnew.end()) << " min dist: " << *std::min_element(distnew.begin(), distnew.end()) << "\n");

	BOOST_MESSAGE(
			"abweichung dlipid new: " << deltaLipidThickness(distnew, d_lipid) << "\n");
	BOOST_MESSAGE(
			"abweichung dlipid oldgrid: " << deltaLipidThickness(distold, d_lipid) << "\n");

	// check if average abbreviation of lipid thickness is below given threshold
	BOOST_CHECK_SMALL(deltaLipidThickness(distold, d_lipid), 10E-6);
}
BOOST_AUTO_TEST_SUITE_END()
