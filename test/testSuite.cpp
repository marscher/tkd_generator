/*
 * testSuite.cpp
 *
 *  Created on: 21.03.2012
 *      Author: marscher
 */
#ifdef BOOST_TEST_DYN_LINK
	#define BOOST_TEST_MAIN
	#include <boost/test/unit_test.hpp>
#else
	#include <boost/test/included/unit_test.hpp>
#endif

#include "common/util/parameter_parsing.h"
#include "validateSwellResults.h"
#include "ug.h"

using namespace boost::unit_test;

struct ug_env {
	ug_env() {
		int argc = framework::master_test_suite().argc;
		char** argv = framework::master_test_suite().argv;
		ug::UGInit(&argc, &argv);
	}
	~ug_env() {
		ug::UGFinalize();
	}
};

BOOST_GLOBAL_FIXTURE(ug_env);

test_suite*
init_unit_test_suite(int argc, char* argv[]) {
	framework::master_test_suite().p_name.set("tkd-Testsuite");

	const char* datafile = NULL;

	if(ug::ParamToString(&datafile, "--datafile", argc, argv)) {
		framework::master_test_suite().add(initValidateResultsTS(datafile));
	}

	return 0;
}
