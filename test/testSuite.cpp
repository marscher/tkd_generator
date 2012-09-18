/*
 * testSuite.cpp
 *
 *  Created on: 21.03.2012
 *      Author: marscher
 */
#include <boost/test/included/unit_test.hpp>
#include <boost/test/results_reporter.hpp>

#include "common/util/parameter_parsing.h"

#include "validateSwellResults.h"
#include "ug.h"

using namespace boost::unit_test;


test_suite*
init_unit_test_suite(int argc, char* argv[]) {
	ug::UGInit(&argc, &argv);

	framework::master_test_suite().p_name.set("tkd-Testsuite");

	const char* datafile = NULL;

	if(ug::ParamToString(&datafile, "--datafile", argc, argv)) {
		framework::master_test_suite().add(initValidateResultsTS(datafile));
	}

	return 0;
}
