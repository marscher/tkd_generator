/*
 * testVolumeOrientation.cpp
 *
 *  Created on: 22.03.2013
 *      Author: marscher
 */
#include <boost/test/unit_test.hpp>
#include "../domain_generator.h"
#include <lib_grid/lib_grid.h>
#include <lib_disc/common/geometry_util.h>
#include <map>
#include <vector>
#include <sstream>
using namespace std;
using namespace ug;

struct gridFixture {
	gridFixture() : grid(),
					sh(grid),
					domGen(grid, sh),
					aaPos(domGen.getVertexAttachmentAccessor())
	{
		grid.attach_to_vertices(aPosition);
	}

	Grid grid;
	SubsetHandler sh;
	tkd::TKDDomainGenerator domGen;
	Grid::VertexAttachmentAccessor<APosition>& aaPos;
};

void performSubsetChecks(Grid&g, SubsetHandler& sh,
		Grid::VertexAttachmentAccessor<APosition>& aaPos,
		const string& sa, const string& sb) {
	int si_a, si_b;
	si_a = sh.get_subset_index(sa.c_str());
	si_b = sh.get_subset_index(sb.c_str());
	BOOST_REQUIRE_MESSAGE(si_a != -1,
			"subset index for subset named " << sa<< " not found");
	BOOST_REQUIRE_MESSAGE(si_b != -1,
			"subset index for subset named " << sb<< " not found");

	GeometricObjectCollection goc1 = sh.get_geometric_objects_in_subset(si_a);
	GeometricObjectCollection goc2 = sh.get_geometric_objects_in_subset(si_b);

	BOOST_CHECK(goc1.num<Face>() == goc2.num<Face>());
	vector3 n1, n2;
	CalculateNormal(n1, *goc1.faces_begin(0), aaPos);
	CalculateNormal(n2, *goc2.faces_begin(0), aaPos);

	// check subsets are anti parallel
	vector3 sum;
	VecAdd(sum, n1, n2);
	BOOST_CHECK_SMALL(sum.x, 1e-6);
	BOOST_CHECK_SMALL(sum.y, 1e-6);
	BOOST_CHECK_SMALL(sum.z, 1e-6);
}

BOOST_FIXTURE_TEST_SUITE(stuff, gridFixture)

BOOST_AUTO_TEST_CASE(test_subset_parallelity) {
	domGen.createSCDomain(10, 30, 24, 5);
	stringstream ss;

	// hex_n_a|b
	// 6 hexas, but two handled per iteration
	for(uint i = 0; i < 3; i++) {
		ss.str("");
		ss << "hex" << i;

		string subset_a = ss.str().append("a");
		string subset_b = ss.str().append("b");
		performSubsetChecks(grid, sh, aaPos, subset_a, subset_b);
	}

	// hex_n_a|b
	// 8 quads (- 2 for top/bottom = 6), but two handled per iteration
	for(uint i = 0; i < 3; i++) {
		ss.str("");
		ss << "quad" << i;

		string subset_a = ss.str().append("a");
		string subset_b = ss.str().append("b");
		performSubsetChecks(grid, sh, aaPos, subset_a, subset_b);
	}

	// perform checks for top and bottom
	performSubsetChecks(grid, sh, aaPos, "top", "bottom");
}

//BOOST_AUTO_TEST_CASE(test_vol_orientation) {
//	domGen.createSCDomain(10, 30, 24, 5);
//	map<Volume*, vector3> vol_normal_map;
//	multimap<Volume*, vector3>::iterator iter;
//	vector<VertexBase*> vertices;

	// iterate over all volumes and check, that normal of base area of each volume is
	// parallel to z axis
//	UG_LOG("fixed orientation: "
//			<< FixOrientation(grid, grid.begin<Volume>(), grid.end<Volume>(),
//					domGen.getVertexAttachmentAccessor()));
//	for(VolumeIterator iter = grid.begin<Volume>(); iter!= grid.end<Volume>(); ++iter) {
//		Volume* v = *iter;
//
//		if(not CheckOrientation(v,domGen.getVertexAttachmentAccessor())) {
//			sh.assign_subset(v, si);
////			FixOrientation()
//		}
//
////		vector3 normal;
////		VolumeVertices::ConstVertexArray arr = v->vertices();
////		for(uint i = 0; i < v->num_vertices(); ++i) {
////			ElementNormal(v->reference_object_id(), normal, arr);
////		}
////		for(uint i = 0; i < v->num_faces(); i++)  {
////			FaceDescriptor fd = v->face_desc(i);
////			CalculateNormal(normal, &fd, domGen.getVertexAttachmentAccessor());
////			vol_normal_map[v] = normal;
////		}
////		FaceDescriptor fd = v->face_desc(0);
////		CalculateNormal(normal, &fd, domGen.getVertexAttachmentAccessor());
////		vol_normal_map[v] = normal;
//	}
////	vector3 zAxis(0,0,1);
////	int si = 7;
////	for(iter = vol_normal_map.begin(); iter != vol_normal_map.end(); ++iter) {
//////		UG_LOG(iter->second << endl);
////		double dot = VecDot(zAxis, iter->second);
////		double deg = rad_to_deg(dot);
////		if(fabs(deg)< 1e-9) {
////			UG_LOG("deg: " << deg << endl);
////			sh.assign_subset(iter->first, si++);
////		}
////	}
//	SaveGridToFile(grid, sh, "/tmp/foo.ugx");
//}

BOOST_AUTO_TEST_SUITE_END();
