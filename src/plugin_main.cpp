/*
 * plugin_main.cpp
 *
 *  Created on: 14.03.2012
 *      Author: marscher
 */

#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "domain_generator.h"
#include "geometry_generator.h"

using namespace ug::bridge;

namespace ug {
namespace tkd {

/** 
 *  \defgroup tkd_generator TKDGenerator
 *  \ingroup plugins_experimental
 *  \brief Generator for Tetrakaidekahedra (TKD).
 *  \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDom>
static void Domain(Registry& reg, std::string grp)
{
	typedef TKDGeometryGenerator geomGen;
	typedef TKDDomainGenerator domGen;

	// register tkd domain generator class
	reg.add_class_<domGen>("TKDDomainGenerator", grp, "Domain (grid) generator for tkds.")
		.add_constructor<void (*)(Grid&, ISubsetHandler&)>(
			"Grid to fill with TKD#SubsetHandler to use")

		.add_constructor<void (*)(Grid&, ISubsetHandler&, bool)>(
			"Grid to fill with TKD"
			"#SubsetHandler to use"
			"#whether a SC domain should be created"
			"#distinct subset for every outer hexagon and quad")

		.add_constructor<void (*)(ug::Domain<3>&)>(
				"fill grid and subset handler of given domain")

		.add_method("setIsSCDomain", &domGen::setSCDomain,
			"switch whether a stratum corneum domain or a simple tkd domain"
			"will be created on createDomain()")

		.add_method("setGridObject", &domGen::setGridObject,
			"sets Grid and SubsetHandler to use.")

	// register createSCDomain
		.add_method(
			"createSCDomain", &domGen::createSCDomain,
			"fills your given grid and SubsetHandler with a stratum corneum domain "
			" with given parameters a, h, w and d_lipid",
			"a#w#h#d_lipid#rows#cols#layers")
	// register createSimpleTKDDomain
		.add_method(
			"createSimpleTKDDomain", &domGen::createSimpleTKDDomain,
			"fills your given grid and SubsetHandler with tkds "
			"with given parameters a, h, w and d_lipid",
			"a#w#h#d_lipid#rows#cols#layers")

	// register TKDGeometryGenerator& getGeometryGenerator() const
		.add_method("getGeometryGenerator",	&domGen::getGeometryGenerator,
			"Geometry generator used to build coordinate informations. \
			 Also calculates volume and surface of single tkd");

	///////////////////////////////////////////////////////////////////////////
	// TKD geometry generator class
	reg.add_class_<geomGen>("TKDGeometryGenerator", grp,
			"Geometry generator class for tkds")
	// ctor
		.add_constructor<void (*)(number, number, number, number)>(
			"height#baseEdgeLength#diameter#d_lipid")
	// register createGeometry()
		.add_method("createGeometry", &geomGen::createGeometry)

	// add volume methods
		.add_method("GetVolume",
			(number (geomGen::*)(int) const) (&geomGen::getVolume),
			"Volume of given subset (1 element)")

	// register surface methods
		.add_method("GetSurface",
			(number (geomGen::*)(int) const) (&geomGen::getSurface),
			"Surface of given subset (1 element)");

	///////////////////////////////////////////////////////////////////////////
	// add static member function of TKDGeometryGenerator as global function
	reg.add_function("GetVolume",
			(number (*)(number, number, number)) (&geomGen::getVolume), grp,
			"calculate Volume for given geometrical parameters");


	// add static member function of TKDGeometryGenerator as global function
	reg.add_function("GetSurface",
			(number (*)(number, number, number)) (&geomGen::getSurface), grp,
			"calculate Surface for given geometrical parameters");

	// register helper functions for subset indices as global functions
	reg.add_function("CorneocyteIndex",
			(int (*)(void)) (&domGen::getCorneocyteIndex), grp,
				"Get subset index for corneocytes");

	reg.add_function("LipidIndex",
			(int (*)(void)) &domGen::getLipidIndex, grp,
				"Get subset index for lipid matrix");
}

}; // end Functionality

// end group tkd_generator
/// \}

} // end namespace tkd

// register tkd generator functions for usage in ug_script
extern "C" void
InitUGPlugin_TKDGenerator(Registry* reg, std::string grp) {
	grp.append("TKDGenerator/");
	RegisterDomain3dDependent<tkd::Functionality>(*reg, grp);
}

} // end namespace ug
