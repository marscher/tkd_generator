/*
 * generator.h
 *
 *  Created on: 08.12.2011
 *      Author: marscher
 */

#ifndef GEOM_GENERATOR_H_
#define GEOM_GENERATOR_H_

#include "common_typedefs.h"
#include "util/rotation_matrix.h"
#include "util/vecComparator.h"
#include <boost/bimap.hpp>

namespace ug {
namespace tkd {

/// \addtogroup tkd_generator
/// \{

/**
 * generates domain decomposition of tetrakaidecahedron.
 * After calling createDomain() you can get the positions and indices
 * of the created geometries to create a Grid object
 */
class TKDGeometryGenerator {
public:

	/**
	 * map indices and coordinates bidirectional (important for domain generator)
	 */
	typedef boost::bimap<int, vector3> CoordIndexMap;

	void createGeometry();

	/**
	 * creates an empty geometry generator. WARNING: this is not yet usable!
	 * use setBaseEdgeLength(), setDiameter(), setHeight() and setLipidThickness()
	 */
	TKDGeometryGenerator();

	/**
	 * @param height
	 * @param baseEdgeLength
	 * @param diameter
	 * @param d_lipid
	 */
	TKDGeometryGenerator(number a, number w, number height, bool createLipid = true,
			number d_lipid = 1);

	/**
	 * should lipid geometry be created?
	 */
	bool createLipid() const;

	void setCreateLipid(bool);

	/**
	 * gets volume of given subset
	 */
	number getVolume(int subset = LIPID) const;

	/**
	 * calculates volume for given geometrical parameters
	 */
	static number getVolume(number a, number s, number h);

	/**
	 * gets surface for given subset
	 */
	number getSurface(int subset = LIPID) const;

	/**
	 * calculates surface for given geometrical parameters
	 */
	static number getSurface(number a, number s, number h);

	CoordIndexMap& getPositions();

	const IndexArray& getIndices() const {
		return indsOut;
	}

	/**
	 * get amount of volumes created. Depends on boolean flag createLipid
	 */
	uint getNumberVolumes() const;

	number getHeight() const;
	number getOverlap() const;

	// setter
	void setHeight(number height);
	void setBaseEdgeLength(number a);
	void setDiameter(number w);
	void setLipidThickness(number d);
	void setGeometricParams(number a, number w, number h, number d_lipid);

protected:
	/**
	 * number of vertices of used solid figures
	 */
	enum SolidFigures {
		Tetrahedron = 4, Pyramid = 5, Prism = 6, Hexahedron = 8
	};

	/**
	 * number of vertices of used solid figures
	 */
	bool b_createLipid;

	// height of one third of tkd : h = 1/3 * h_tkd
	number h_corneocyte;
	// base edge length of central hexahedron
	number a_corneocyte;
	// diameter
	number w_corneocyte;
	number s_corneocyte;
	// thickness of lipid matrix
	number d_lipid;
	// base edge length of lipid matrix
	number a_lipid;
	// height of 1/3 of lipid matrix
	number h_lipid;
	// quantity s, overlap of two aligned tkds with lipid matrix
	number s_lipid;
	number w_lipid;

	/**
	 * stores coordinates of points
	 */
	CoordIndexMap m_mCoordIndexMapping;

	/**
	 * stores indices for points
	 */
	IndexArray indsOut;

	// matrix used to rotate all geometric objects around z axis
	RotationMatrix R;

	/**
	 * count how much unique vertices have been stored in output map
	 */
	uint m_num_verts_created;

	/**
	 * inits base geometric parameters
	 */
	void initGeometricParams();

	void initLipidGeometric();

	void updateOverlap(int subset);

	void setLipidBaseEdgeLength();
	void setLipidHeight();
	void setLipidDiameter();
	void setLipidParameters();

	void createCorneocyte(const vector3&);
	void createLipidMatrix(const vector3&, const number rotationOffset = 0,
			bool bottom = false);

	/**
	 *	creates the upper part of tkd (symmetric to bottom part!)
	 */
	void createCorneocyteTop(const vector3& offset,
			const number rotationOffset = 0, bool bottom = false);

	/**
	 * creates middle part with given offset
	 */
	void createCorneocyteMiddle(const vector3& offset);

	/**
	 * pushes posIn into global posOut reference and creates
	 * indices for each vertex which are pushed into global indsOut reference
	 */
	void createGeometricObject(const CoordsArray& posIn);

	void flipOrientationPrism(CoordsArray& prismPos) const;
	void flipOrientationPyramid(CoordsArray& pyramidPos) const;
	void flipOrientationTetrahedron(CoordsArray& tetrahedronPos) const;
	void flipOrientationHexahedron(CoordsArray& hexaPos) const;

	///// segments of top and bottom
	CoordsArray topInner;
	CoordsArray topOuterPrism;

	CoordsArray topOuter_Pr_rightTetrahedron;
	CoordsArray topOuter_Pr_leftTetrahedron;
	CoordsArray topOuter_Pr2T_prism;

	///// segments of middle part
	// outer prism
	CoordsArray middleOuterPrism1;
	// inner prism
	CoordsArray middleOuterPrism2;

	// below obenAussenPr_rightTetrahedron
	CoordsArray middleOuterH2Pr_tetrahedron;
	CoordsArray middleOuterH2Pr_pyramid;

	CoordsArray middleOuter2PrH_tetrahedron;
	CoordsArray middleOuter2PrH_pyramid;

	// below obenAussenPrism
	CoordsArray middleOuterHexahedron;

	// lipid matrix
	/// lipid hexagon
	CoordsArray topInnerPrismL;
	CoordsArray topOuterPrismL;
	CoordsArray topOuter_leftPrismL;
	CoordsArray topOuter_rightPrismL;
	CoordsArray upperHexahedronL;
	CoordsArray bottomOuterHexahedronL;
	CoordsArray bottomLeftPrismL;
	CoordsArray bottomRightPrismL;

	// lipid side quad
	CoordsArray sideQuad;

	const static vector3 origin;
};

// end group tkd_generator
/// \}

} //end of namespace tkd
}// end of namespace ug
#endif /* GEOM_GENERATOR_H_ */
