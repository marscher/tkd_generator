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
#include "domain_generator.h"

#include <utility>
#include <cmath>

using namespace boost::unit_test;
using namespace ug;
using namespace ug::tkd;

BOOST_AUTO_TEST_SUITE(lipid)

BOOST_AUTO_TEST_CASE(compareThicknessWithOld) {
	double a = 10;
	double h = sqrt(6) * a;
	double w = 3 * a;
	double d_lipid = 8;

	// init libgrid objects
	Grid gridNew(GRIDOPT_STANDARD_INTERCONNECTION);
	Grid gridOld(GRIDOPT_STANDARD_INTERCONNECTION);
	SubsetHandler shold(gridOld);
	SubsetHandler shnew(gridNew);

	// create grid with new tkd generator
	gridNew.attach_to_vertices(aPosition);
	gridOld.attach_to_vertices(aPosition);

	TKDDomainGenerator gen(gridNew, shnew);
	gen.createSCDomain(a, w, h, d_lipid);

	// grid created with same params with old tkdmodeller
	std::string oldfile = "nucleus3d_10-24.495-30-8.000_1x1x1.ugx";

	bool gridLoaded = LoadGridFromFile(gridOld, shold, oldfile.c_str());
	BOOST_REQUIRE_MESSAGE(gridLoaded,
			"old grid file " << oldfile << " could not be loaded.");

	vector<number> vols_old = meassureLipidVolume(gridOld, shold, d_lipid);
	vector<number> vols_new = meassureLipidVolume(gridNew, shnew, d_lipid);

	// check if average abbreviation of lipid thickness is below given threshold
	double sum_volume_old = std::accumulate(vols_old.begin(), vols_old.end(), 0);
	double sum_volume_new = std::accumulate(vols_new.begin(), vols_new.end(), 0);

	BOOST_REQUIRE_CLOSE(sum_volume_new, sum_volume_old, 10E-16);
}
BOOST_AUTO_TEST_SUITE_END()
