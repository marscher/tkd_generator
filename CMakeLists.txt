################################################################################
# TKD_GENERATOR
#
# created by Martin Scherer
# 
# This plugin is used to create a domain consisting of 
# overlapping tetrakaidecahedrons
################################################################################

# Name of your plugin, path to your plugin (relative to ug's root path) and sources.
set(pluginName	TKDGenerator)
set(SOURCES src/plugin_main.cpp
            src/domain_generator.cpp
            src/geometry_generator.cpp
            src/coordinates.cpp
            src/util/rotation_matrix.cpp
            src/util/vecComparator.cpp
            src/util/geometric_helper.cpp)
			
# sources test
set(SOURCES_TEST test/testSuite.cpp 
                 test/testLipidDistance.cpp
                 test/testVolume.cpp
                 test/testVolumeOrientation.cpp
                 test/validateSwellResults.cpp
                 src/util/volume_calculation.cpp)
				 
option(TKDTestsuite "Build TKDTestsuite" OFF)
message(STATUS "Info: TKDTestsuite: ${TKDTestsuite} (options are: ON, OFF)")

cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})
# include the definitions and dependencies for ug-plugins.
include("$ENV{UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake")
# needed for initialization of ug parallel env in boost testing framework
include_directories("$ENV{UG_ROOT_PATH}/unit_tests/src" src)

if(NOT "${DIM}" STREQUAL "ALL")
	list(FIND DIM "3" 3d_enabled)
	IF(3d_enabled EQUAL -1)
		message(WARNING "${pluginName} plugin requires DIM=3, DIM=ALL, DIM=x;3")
	ENDIF()
endif()


if(TKDTestsuite)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)

	find_package(Boost COMPONENTS unit_test_framework)
	if(Boost_UNIT_TEST_FRAMEWORK_LIBRARY)
		add_definitions(-DBOOST_TEST_DYN_LINK)
	endif()

	add_executable(TKDTestsuite ${SOURCES_TEST})
endif(TKDTestsuite)

################################################################################
# The code below doesn't have to be changed (usually)
################################################################################

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${pluginPath} ${SOURCES})
	if(TKDTestsuite)
		target_link_libraries (TKDTestsuite ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ug4_s)
	endif(TKDTestsuite)
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
	if(TKDTestsuite)
		target_link_libraries (TKDTestsuite TKDGenerator ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ug4)
	endif(TKDTestsuite)
endif(buildEmbeddedPlugins)
