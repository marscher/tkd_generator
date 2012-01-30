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

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")
# needed for initialization of ug parallel env in boost testing framework
include_directories("${UG_ROOT_PATH}/unit_tests/src")

# sources plugin
set(SOURCES	tetrakaidecahedron_generator.cpp
			rotation_matrix.cpp
			geometric_helper.cpp
			generator.cpp
			util/volume_calculation.cpp)

# sources test
set(SOURCES_TEST test/testLipidDistance.cpp
				 test/testVolume.cpp)

# create a shared library from the sources and link it against ug4.
add_library(TKDGenerator SHARED ${SOURCES})
# TODO: dont write exe to plugin out dir
add_executable(TKDTestsuite ${SOURCES_TEST})

if(STATIC)
	target_link_libraries (TKDGenerator ug4_s)
	target_link_libraries (TKDTestsuite TKDGenerator ug4_s)
else(STATIC)
	target_link_libraries (TKDGenerator ug4)
	target_link_libraries (TKDTestsuite TKDGenerator ug4)
endif(STATIC)