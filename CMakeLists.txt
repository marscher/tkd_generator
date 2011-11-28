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

# sources
set(SOURCES	tetrakaidecahedron_generator.cpp)


# create a shared library from the sources and link it against ug4.
add_library(TKDGenerator SHARED ${SOURCES})

if(STATIC)
	target_link_libraries (SamplePlugin ug4_s)
else(STATIC)
	target_link_libraries (SamplePlugin ug4)
endif(STATIC)