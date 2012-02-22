#ifndef __coordinates_h__
#define __coordinates_h__

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

namespace tkdGenerator {
using std::vector;

inline vector<vector3> CalculateLipidCoords(number a_corneo, number high,
		number width, number d_lipid, const vector3& offset) {
	number b = sqrt(3);

	// lipid coords
	vector<vector3> l(38);
	// corneocyte coords
	vector<vector3> c(76);

	unsigned int i;

	number a1 = sqrt(
			1. / 9. * high * high
					+ 1. / 3. * (width - 2. * a_corneo) * (width - 2. * a_corneo));

	number alpha = acos((width - 2.0 * a_corneo) / (2.0 * a1));
	number beta = 90. / 180. * M_PI + acos(1. / 3. * high / (a1 * sin(alpha)));
	number gamma = acos(1. / 3. * high / a1) + 90. / 180. * M_PI;

	number high_l = high + d_lipid;

	number dis = d_lipid / 2.0;

	number m1 = dis / tan(beta / 2.0);
	number m2 = dis / tan(gamma / 2.0);
	number a_lipid = (b * a_corneo + m1 + m2) / b;

	number a2 = 1. / 3. * high_l / sin(M_PI - gamma);

	number bc = a1 * cos(alpha);

	// x offset von mitte, welches element???
	number bc1 = a2 * cos(alpha);

	number ab = a1 * cos(alpha) * 2.0 / b;
	number ab1 = a2 * cos(alpha) * 2.0 / b;

	number ac = a1 * cos(alpha) / b;
	number ac1 = a2 * cos(alpha) / b;

	number ec = ac * b / 2.0;
	number ec1 = ac1 * b / 2.0;
	number af = ac * 0.5;
	number af1 = ac1 * 0.5;

	l[0] = vector3(-a_lipid / 2.0, -a_lipid * b / 2.0, dis);
	l[1] = vector3(-l[0].x, l[0].y, l[0].z);
	l[2] = vector3(a_lipid, 0, l[0].z);
	l[3] = vector3(l[1].x, -l[1].y, l[0].z);
	l[4] = vector3(l[0].x, -l[0].y, l[0].z);
	l[5] = vector3(-l[2].x, l[2].y, l[0].z);
	l[6] = vector3(0, 0, l[0].z);

	c[0] = vector3(-a_corneo / 2.0, -a_corneo * b / 2.0, 0);
	c[1] = vector3(-c[0].x, c[0].y, c[0].z);
	c[2] = vector3(a_corneo, 0, c[0].z);
	c[3] = vector3(c[1].x, -c[1].y, c[0].z);
	c[4] = vector3(c[0].x, -c[0].y, c[0].z);
	c[5] = vector3(-c[2].x, c[2].y, c[0].z);
	c[6] = vector3(0, 0, 0);

	for (i = 0; i < 7; i++) {
		l[i + 7] = c[i];
	}

	c[7] = vector3(c[0].x, c[0].y - ab, c[0].z - 1. / 3 * high);

	l[14] = vector3(l[0].x, l[0].y - ab1, l[0].z - 1.0 / 3.0 * high_l);
	l[15] = vector3(l[1].x, l[14].y, l[14].z);
	l[16] = vector3(l[1].x + ec1, l[1].y - af1, l[14].z);
	l[17] = vector3(l[2].x + ec1, l[2].y - af1, l[14].z);
	l[18] = vector3(l[2].x + bc1, l[2].y + ac1, l[14].z);
	l[19] = vector3(l[3].x + bc1, l[3].y + ac1, l[14].z);
	l[20] = vector3(l[3].x, l[19].y, l[14].z);
	l[21] = vector3(l[4].x, l[19].y, l[14].z);
	l[22] = vector3(-l[19].x, l[19].y, l[14].z);
	l[23] = vector3(-l[18].x, l[18].y, l[14].z);
	l[24] = vector3(-l[17].x, l[17].y, l[14].z);
	l[25] = vector3(-l[16].x, l[16].y, l[14].z);

	c[8] = vector3(c[1].x, c[7].y, c[0].z - 1.0 / 3.0 * high);
	c[9] = vector3(c[1].x + ec, c[1].y - af, c[8].z);
	c[10] = vector3(c[2].x + ec, c[2].y - af, c[9].z);
	c[11] = vector3(c[2].x + bc, c[2].y + ac, c[10].z);
	c[12] = vector3(c[3].x + bc, c[3].y + ac, c[11].z);
	c[13] = vector3(c[3].x, c[12].y, c[12].z);
	c[14] = vector3(c[4].x, c[13].y, c[13].z);
	c[15] = vector3(-c[12].x, c[14].y, c[14].z);
	c[16] = vector3(-c[11].x, c[11].y, c[15].z);
	c[17] = vector3(-c[10].x, c[10].y, c[16].z);
	c[18] = vector3(-c[9].x, c[9].y, c[17].z);

	for (i = 7; i < 19; i++) {
		l[19 + i] = c[i];
	}

	c[31] = vector3(c[0].x, c[0].y, c[0].z - high);
	c[32] = vector3(c[1].x, c[1].y, c[31].z);
	c[33] = vector3(c[2].x, c[2].y, c[31].z);
	c[34] = vector3(c[3].x, c[3].y, c[31].z);
	c[35] = vector3(c[4].x, c[4].y, c[31].z);
	c[36] = vector3(c[5].x, c[5].y, c[31].z);
	c[37] = vector3(c[6].x, c[6].y, c[31].z);

	l[69] = vector3(l[0].x, l[0].y, l[0].z - high_l);
	l[70] = vector3(l[1].x, l[1].y, l[69].z);
	l[71] = vector3(l[2].x, l[2].y, l[69].z);
	l[72] = vector3(l[3].x, l[3].y, l[69].z);
	l[73] = vector3(l[4].x, l[4].y, l[69].z);
	l[74] = vector3(l[5].x, l[5].y, l[69].z);
	l[75] = vector3(l[6].x, l[6].y, l[69].z);

	c[34] = vector3(c[3].x, c[3].y, c[31].z);
	c[35] = vector3(c[4].x, c[4].y, c[31].z);
	c[36] = vector3(c[5].x, c[5].y, c[31].z);
	c[37] = vector3(c[6].x, c[6].y, c[31].z);

	for (i = 31; i < 38; i++) {
		l[31 + i] = c[i];
	}

	c[19] = vector3(c[31].x - bc, c[31].y - ac, c[0].z - 2.0 / 3.0 * high);
	c[20] = vector3(c[31].x, c[19].y, c[19].z);
	c[21] = vector3(c[32].x, c[20].y, c[20].z);
	c[22] = vector3(-c[19].x, c[21].y, c[20].z);
	c[23] = vector3(c[33].x + bc, c[33].y - ac, c[22].z);

	c[24] = vector3(c[33].x + ec, c[33].y + af, c[22].z);
	c[25] = vector3(c[34].x + ec, c[34].y + af, c[24].z);
	c[26] = vector3(c[34].x, c[34].y + ab, c[25].z);
	c[27] = vector3(c[35].x, c[26].y, c[26].z);
	c[28] = vector3(-c[25].x, c[25].y, c[27].z);
	c[29] = vector3(-c[24].x, c[24].y, c[28].z);
	c[30] = vector3(-c[23].x, c[23].y, c[29].z);

	l[50] = vector3(l[69].x - bc1, l[69].y - ac1, l[0].z - 2.0 / 3.0 * high_l);
	l[51] = vector3(l[69].x, l[50].y, l[50].z);
	l[52] = vector3(l[70].x, l[50].y, l[50].z);
	l[53] = vector3(-l[50].x, l[50].y, l[50].z);
	l[54] = vector3(l[71].x + bc1, l[71].y - ac1, l[50].z);
	l[55] = vector3(l[71].x + ec1, l[71].y + af1, l[50].z);
	l[56] = vector3(l[72].x + ec1, l[72].y + af1, l[50].z);
	l[57] = vector3(l[72].x, l[72].y + ab1, l[50].z);
	l[58] = vector3(l[73].x, l[57].y, l[50].z);
	l[59] = vector3(-l[56].x, l[56].y, l[50].z);
	l[60] = vector3(-l[55].x, l[55].y, l[50].z);
	l[61] = vector3(-l[54].x, l[54].y, l[50].z);

	for (i = 19; i < 31; i++) {
		l[19 + i] = c[i];
	}

	for (i = 0; i < l.size(); i++) {
		l[i].x += offset.x;
		l[i].y += offset.y;
		l[i].z += offset.z;
	}

	return l;
}
}	// end of namespace tkdGenerator
#endif
