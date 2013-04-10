################################################################################
# TKD_GENERATOR
#
# created by Martin Scherer
# 
# This plugin is used to create a domain consisting of 
# overlapping tetrakaidecahedrons
################################################################################

cmake_minimum_required(VERSION 2.6)

project(TKD_GENERATOR)

if(NOT "${DIM}" STREQUAL "ALL")
	list(FIND DIM "3" 3d_enabled)
	IF(3d_enabled EQUAL -1)
		message(WARNING "TKDGenerator plugin requires DIM=3, DIM=ALL, DIM=x;3")
	ENDIF()
endif()

option(TKDTestsuite "Build TKDTestsuite" OFF)
message(STATUS "Info: TKDTestsuite:	"${TKDTestsuite}" (options are: ON, OFF)")

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")
# needed for initialization of ug parallel env in boost testing framework
include_directories("${UG_ROOT_PATH}/unit_tests/src")

# sources plugin
set(SOURCES plugin_main.cpp
			domain_generator.cpp
			geometry_generator.cpp
			util/rotation_matrix.cpp
			util/vecComparator.cpp
			util/geometric_helper.cpp)

# sources test
set(SOURCES_TEST test/testSuite.cpp
				 test/testLipidDistance.cpp
				 test/testVolume.cpp
				 test/testVolumeOrientation.cpp
				 test/validateSwellResults.cpp
				 util/volume_calculation.cpp)

if(TKDTestsuite)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)

	find_package(Boost COMPONENTS unit_test_framework)
	if(Boost_UNIT_TEST_FRAMEWORK_LIBRARY)
		add_definitions(-DBOOST_TEST_DYN_LINK)
	endif()

	add_executable(TKDTestsuite ${SOURCES_TEST})
endif(TKDTestsuite)

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/tkd_generator ${SOURCES})
	if(TKDTestsuite)
		target_link_libraries (TKDTestsuite
							${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
							ug4_s)
	endif(TKDTestsuite)
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(TKDGenerator SHARED ${SOURCES})
	target_link_libraries (TKDGenerator ug4)
	if(TKDTestsuite)
		target_link_libraries (TKDTestsuite
								TKDGenerator
								${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
								ug4)
	endif(TKDTestsuite)
endif(buildEmbeddedPlugins)
