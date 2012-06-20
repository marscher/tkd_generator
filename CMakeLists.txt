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

option(TKDTestsuite "Build TKDTestsuite" OFF)
message(STATUS "Info: TKDTestsuite:      "${TKDTestsuite}" (options are: ON, OFF)")

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")
# needed for initialization of ug parallel env in boost testing framework
include_directories("${UG_ROOT_PATH}/unit_tests/src")

# sources plugin
set(SOURCES	plugin_main.cpp
            domain_generator.cpp
            geometry_generator.cpp
			util/rotation_matrix.cpp
			util/geometric_helper.cpp
			util/volume_calculation.cpp
			util/vecComparator.cpp)

# sources test
set(SOURCES_TEST test/testSuite.cpp
                 test/testLipidDistance.cpp
				 test/testVolume.cpp
				 test/validateSwellResults.cpp
				 util/VTKOutputGrid.cpp)

if(TKDTestsuite)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)
	add_executable(TKDTestsuite ${SOURCES_TEST})
endif(TKDTestsuite)

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/tkd_generator ${SOURCES})
	if(TKDTestsuite)
		target_link_libraries (TKDTestsuite ug4_s)
	endif(TKDTestsuite)
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(TKDGenerator SHARED ${SOURCES})
	target_link_libraries (TKDGenerator ug4)
	if(TKDTestsuite)
		target_link_libraries (TKDTestsuite TKDGenerator ug4)
	endif(TKDTestsuite)
endif(buildEmbeddedPlugins)