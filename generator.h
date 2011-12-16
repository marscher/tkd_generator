/*
 * generator.h
 *
 *  Created on: 08.12.2011
 *      Author: marscher
 */

#ifndef GENERATOR_H_
#define GENERATOR_H_

#include "geometric_helper.h"
#include "common_typedefs.h"
#include "rotation_matrix.h"

namespace tkdGenerator {

/**
 * number of vertices of used solid figures
 */
enum SolidFigures {
	Tetrahedron = 4, Pyramid = 5, Prism = 6, Hexahedron = 8
};

const static v origin(0, 0, 0);

class Generator {
public:
	void createTKD();
	void createTKD(const v& origin);

	/**
	 * @param height
	 * @param baseEdgeLength
	 * @param diameter
	 * @param posOut
	 * @param indsOut
	 */
	Generator(number height, number baseEdgeLength, number diameter,
			CoordsArray& posOut, IndexArray& indsOut) :
			h(height / 3), a(baseEdgeLength), w(diameter), posOut(posOut), indsOut(
					indsOut), R(0) {
		init();
	}

	number getVolume() const;

	number getSurface() const;

protected:
	// height of one third of tkd : h = 1/3 * h_tkd
	number h;
	// base edge length of central hexahedron
	number a;
	// diameter
	number w;
	// quantity s, overlap of two aligned tkd's
	number s;
	// height of base triangle of top inner prism
	number g;
	// height of base triangle of tetrahedron of ObenAussenPr2T
	number b;

	number centerOfMass;

	/**
	 * stores a reference to CoordsArray owned by creator of this instance
	 */
	CoordsArray& posOut;
	/**
	 * stores a reference to IndexArray owned by creator of this instance
	 */
	IndexArray& indsOut;

	// matrix used to rotate all geometric objects
	RotationMatrix R;

	// global index which is incremented for all vertices generated by this instance
	size_t index;

	/**
	 * inits base geometric parameters
	 */
	void init();

	/**
	 *	creates the upper part of tkd (symmetric to bottom part!)
	 */
	void createTop(const v& offset, const number rotationOffset = 0);

	/**
	 * creates middle part with given offset
	 */
	void createMiddle(const v& offset);

	/**
	 * pushes posIn into global posOut reference and creates
	 * indices for each vertex which are pushed into global indsOut reference
	 */
	void createGeometricObject(const CoordsArray & posIn);

	///// segments of top and bottom
	CoordsArray obenInnen;
	// G(Ki -> ObenAussenPr)
	CoordsArray obenAussenPrism;

	CoordsArray obenAussenPr_rightTetrahedron;
	CoordsArray obenAussenPr_leftTetrahedrson;

	CoordsArray obenAussenPr2T_prism;

	///// segments of middle part
	// outer prism
	CoordsArray mitteAussenP1;
	// inner prism
	CoordsArray mitteAussenP2;

	// symetric with mittAussenH2Pr!
	//	CoordsArray mittAussen2PrH;

	// below obenAussenPr_rightTetrahedron
	CoordsArray mitteAussenH2Pr_tetrahedron;
	CoordsArray mitteAussenH2Pr_pyramid;

	CoordsArray mitteAussen2PrH_tetrahedron;
	CoordsArray mitteAussen2PrH_pyramid;

	// below obenAussenPrism
	CoordsArray mitteAussenHexahedron;
};

} //end of namespace
#endif /* GENERATOR_H_ */
