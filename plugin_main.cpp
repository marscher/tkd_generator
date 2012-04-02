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
	domgenC.add_method<
			void (domGen::*)(number, number, number, number, int, int, int)>(
			"CreateDomain", &domGen::createTKDDomain,
			"fills your given grid and SubsetHandler with tkds with given parameters a, h, w and d_lipid",
			"a#w#h#d_lipid#rows#cols#layers");

	// register setSubsetHandlerInfo()
	domgenC.add_method<
			void (domGen::*)(const char*, const char*, const vector4&,
					const vector4&)>("SetSubsetInfo",
			&domGen::setSubsetHandlerInfo,
			"set SubsetHandler informations (name, color) for corneocytes and lipid matrix",
			"corneocyte_name#lipid_name#corneocyte_color#lipid_color");

	// todo register helper functions for subset indices
//	domgenC.add_method("getCorneocyteIndex",
//				&domGen::getCorneocyteIndex,
//				"Get subset index for corneocytes");
//
//	domgenC.add_method("getLipidIndex",
//				&domGen::getLipidIndex,
//				"Get subset index for lipid matrix");

	// register TKDGeometryGenerator& getGeometryGenerator() const
	domgenC.add_method<geomGen& (domGen::*)(void)>("GetGeometryGenerator",
			&domGen::getGeometryGenerator,
			"geometry generator used to build coordinate informations. \
			 Also calcs volume and surface of single tkd");

	// tkd geometry generator class
	bridge::ExportedClass<geomGen>& geomGenC = reg->add_class_<geomGen>(
			"TKDGeometryGenerator", grp, "Geometry generator class for tkds");
	// ctor
	geomGenC.add_constructor<void (*)(number, number, number, number)>(
			"height#baseEdgeLength#diameter#d_lipid");
	// register createGeometry()
	geomGenC.add_method("createGeometry", &geomGen::createGeometry);

	// add overloaded methods in int subset and number a, h, s
	geomGenC.add_method("getVolume",
			(number (geomGen::*)(int) const)(&geomGen::getVolume), "volume of given subset (1 element)");
	geomGenC.add_method("getVolume",
			(number (geomGen::*)(number, number,
					number) const)(&geomGen::getVolume), "volume for given geometrical parameters");

	// register surface methods
	geomGenC.add_method("getSurface",
			(number (geomGen::*)(int) const)(&geomGen::getSurface), "surface of given subset (1 element)");
	geomGenC.add_method("getSurface",
			(number (geomGen::*)(number, number,
					number) const)(&geomGen::getSurface), "surface for given geometrical parameters");
}
