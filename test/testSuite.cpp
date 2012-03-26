/*
 * testSuite.cpp
 *
 *  Created on: 21.03.2012
 *      Author: marscher
 */
#include <boost/test/included/unit_test.hpp>
#include <boost/test/results_reporter.hpp>

#include "common/util/parameter_parsing.h"
#include "fixtures/UGScriptingEnvFixture.h"

#include "validateSwellResults.h"
//#include <iostream>
//using namespace std;

using namespace boost::unit_test;

//BOOST_GLOBAL_FIXTURE(UGScriptingEnvFixture);

test_suite*
init_unit_test_suite(int argc, char* argv[]) {
	framework::master_test_suite().p_name.value = "tkd-Testsuite";

	int i = ug::GetParamIndex("--datafile", argc, argv);
	if (i > 0) {
		const char* datafile = argv[i + 1];
		framework::master_test_suite().add(validateResultsTS(datafile));
	}
	return 0;
}
