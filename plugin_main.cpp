/*
 * plugin_main.cpp
 *
 *  Created on: 14.03.2012
 *      Author: marscher
 */

#include "registry/registry.h"
#include "bridge/bridge.h"
#include "lib_grid/lib_grid.h"
#include "common/types.h"

#include "domain_generator.h"
#include "geometry_generator.h"

using namespace ug;

// register tkd generator functions for usage in ug_script
extern "C" void InitUGPlugin(bridge::Registry* reg, std::string parentGroup) {
	UG_SET_DEBUG_LEVEL(LogAssistant::APP, 0);

	std::string grp(parentGroup);
	grp.append("tkd_generator/");

	typedef tkd::TKDGeometryGenerator geomGen;
	typedef tkd::TKDDomainGenerator domGen;

	// register tkd domain generator class
	bridge::ExportedClass<domGen>& domgenC = reg->add_class_<domGen>(
			"TKDDomainGenerator", grp, "Domain (grid) generator for tkds.");

	domgenC.add_constructor<void (*)(Grid&, SubsetHandler&)>(
			"Grid to fill with TKD#SubsetHandler to use");

	// register createTKDDomain
	domgenC.add_method(
			"createDomain", &domGen::createTKDDomain,
			"fills your given grid and SubsetHandler with tkds with given parameters a, h, w and d_lipid",
			"a#w#h#d_lipid#rows#cols#layers");

	// register setSubsetHandlerInfo()
	domgenC.add_method("setSubsetInfo",
			&domGen::setSubsetHandlerInfo,
			"Set SubsetHandler informations (name, color) for corneocytes and lipid matrix",
			"corneocyte_name#lipid_name#corneocyte_color#lipid_color");

	// register helper functions for subset indices as global functions
	reg->add_function("CorneocyteIndex",
			(int (*)(void)) (&domGen::getCorneocyteIndex),
				"Get subset index for corneocytes");

	reg->add_function("LipidIndex",
			(int (*)(void)) &domGen::getLipidIndex,
				"Get subset index for lipid matrix");

	// register TKDGeometryGenerator& getGeometryGenerator() const
	domgenC.add_method("getGeometryGenerator",
			&domGen::getGeometryGenerator,
			"Geometry generator used to build coordinate informations. \
			 Also calculates volume and surface of single tkd");

	// tkd geometry generator class
	bridge::ExportedClass<geomGen>& geomGenC = reg->add_class_<geomGen>(
			"TKDGeometryGenerator", grp, "Geometry generator class for tkds");
	// ctor
	geomGenC.add_constructor<void (*)(number, number, number, number)>(
			"height#baseEdgeLength#diameter#d_lipid");
	// register createGeometry()
	geomGenC.add_method("createGeometry", &geomGen::createGeometry);

	// add volume methods
	geomGenC.add_method("getVolume",
			(number (geomGen::*)(int) const) (&geomGen::getVolume), "Volume of given subset (1 element)");
	// add static member function of TKDGeometryGenerator as global function
	reg->add_function("GetVolume",
			(number (*)(number, number, number)) (&geomGen::getVolume), "calculate Volume for given geometrical parameters");

	// register surface methods
	geomGenC.add_method("getSurface",
			(number (geomGen::*)(int) const) (&geomGen::getSurface), "Surface of given subset (1 element)");
	// add static member function of TKDGeometryGenerator as global function
	reg->add_function("GetSurface",
			(number (*)(number, number, number)) (&geomGen::getSurface), "calculate Surface for given geometrical parameters");
}
