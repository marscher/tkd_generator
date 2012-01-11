/*
 * testLipidDistance.cpp
 *
 *  Created on: 10.01.2012
 *      Author: marscher
 */
#define BOOST_TEST_MODULE testLipid

#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "common/util/path_provider.h"
#include "common/os_dependent/file_util.h"

#include "fixtures/UGScriptingEnvFixture.h"
#include "testHelper.h"
#include "../generator.h"

using namespace boost::unit_test;
using namespace tkdGenerator;

BOOST_GLOBAL_FIXTURE(UGScriptingEnvFixture)

BOOST_AUTO_TEST_SUITE(lipid)

/**
 * performs checks if lipid thickness (= distance all boundary faces to inner faces) is correct
 */
void checkThickness(number a, number h, number w, number d_lipid) {
//	number a = 10;
//	number h = sqrt(6)*a;
//	number w = 3*a;
//	number d_lipid = 8;

// init libgrid objects
	Grid gridNew(GRIDOPT_STANDARD_INTERCONNECTION);
	SubsetHandler shnew(gridNew);
	vector<number> distnew = meassureLipidThickness(gridNew, shnew, d_lipid);

// create grid with new tkd generator
	TKDGeometryGenerator gen(h, a, w, d_lipid);
	gen.createDomain();
	gridNew.attach_to_vertices(aPosition);
	createGrid(gridNew, shnew, gen.getPositions(), gen.getIndices());

//	BOOST_MESSAGE(
//			"max dist: " << *std::max_element(distnew.begin(), distnew.end())
//			<< " min dist: " << *std::min_element(distnew.begin(), distnew.end()) << "\n");
//
//	BOOST_MESSAGE(
//			"abweichung dlipid new: " << deltaLipidThickness(distnew, d_lipid) << "\n");

// check if average abbreviation of lipid thickness is below given threshold
	BOOST_CHECK_SMALL(deltaLipidThickness(distnew, d_lipid), 10E-6);
}

void foo(int a, int b) {
	BOOST_CHECK( a + b < 20);
}

struct sub_test_suite : public test_suite {

	sub_test_suite() : test_suite("foo") {
		parameters_list.push_back( 1 );
		parameters_list.push_back( 5 );
		parameters_list.push_back( 6 );
		parameters_list.push_back( 7 );
		parameters_list.push_back( 140 );
		number a = 10;
		number h = sqrt(6)*a;
		number w = 3*a;
		number d_lipid = 8;

//		add(
//				BOOST_PARAM_TEST_CASE(
//						boost::callback3<number, number, number>( boost::bind(&foo, _1, h, w, a)),
//						parameters_list.begin(), parameters_list.end()));
//		boost::function2<number,number,number>(boost::bind( &checkThickness(), _1, h, w, d_lipid));
//		add( BOOST_PARAM_TEST_CASE( f,
//						parameters_list.begin(), parameters_list.end() ));
		framework::master_test_suite().add(this);
	}

	std::list<number> parameters_list;
	std::list<number> rangeA;
	std::list<number> rangeH;
};

BOOST_AUTO_TEST_CASE(compareThicknessWithOld) {
	number a = 10;
	number h = sqrt(6)*a;
	number w = 3*a;
	number d_lipid = 8;

	// init libgrid objects
	Grid gridNew(GRIDOPT_STANDARD_INTERCONNECTION);
	Grid gridOld(GRIDOPT_STANDARD_INTERCONNECTION);
	SubsetHandler shold(gridOld);
	SubsetHandler shnew(gridNew);

	// create grid with new tkd generator
	TKDGeometryGenerator gen(h, a, w, d_lipid);
	gen.createDomain();
	gridNew.attach_to_vertices(aPosition);
	createGrid(gridNew, shnew, gen.getPositions(), gen.getIndices());

	// grid created with same params with old tkdmodeller
	string oldfile = PathProvider::get_path(ROOT_PATH)
	+ "/plugins/experimental/tkd_generator/test/nucleus3d_10-24.495-30-8.000_1x1x1.ugx";

	bool gridLoaded = LoadGridFromFile(gridOld, shold, oldfile.c_str());
//	BOOST_REQUIRE_MESSAGE(gridLoaded, "old grid file "
//			<< oldfile << " could not be loaded.");

	vector<number> distold = meassureLipidThickness(gridOld, shold, d_lipid);

	BOOST_MESSAGE("______________________________________\n");
	vector<number> distnew = meassureLipidThickness(gridNew, shnew, d_lipid);

	BOOST_MESSAGE(
			"max dist: " << *std::max_element(distold.begin(), distold.end())
			<< " min dist: " << *std::min_element(distold.begin(), distold.end()) << "\n" );
	BOOST_MESSAGE("______________________________________\n");
	BOOST_MESSAGE(
			"max dist: " << *std::max_element(distnew.begin(), distnew.end())
			<< " min dist: " << *std::min_element(distnew.begin(), distnew.end()) << "\n");

	BOOST_MESSAGE(
			"abweichung dlipid new: " << deltaLipidThickness(distnew, d_lipid) << "\n");
	BOOST_MESSAGE(
			"abweichung dlipid oldgrid: " << deltaLipidThickness(distold, d_lipid) << "\n");

// check if average abbreviation of lipid thickness is below given threshold
	BOOST_CHECK_SMALL(deltaLipidThickness(distold, d_lipid), 10E-6 );
}

BOOST_AUTO_TEST_SUITE_END()
