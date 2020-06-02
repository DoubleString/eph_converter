/*
 * Broadcast.cpp
 *
 *  Created on: 2020年5月11日
 *      Author: dongdong
 */

#include <iostream>
#include <cmath>
#include "../include/Broadcast.h"
#include "../include/Com.h"

using namespace std;

int BroadcastEphUtils::m_brd2xyz(const char* mode, const char* cprn, int wk,
		double sow, double* xsat, double* clk, double* dtmin, double* tgd,
		int* iode) {
	int i, isys, lfind,f_prn;
	double GMEARTH, EARTH_ROTATE;
	double dt, dtc, pos[3], svpos[3], geovel[3], detavel[3];
	double a, xn, xm, ex, e, v0, vs, vc, phi, ccc, sss, du, dr, di, r, u, xi,
			xx, yy, intv;
	double xnode, xpdot, ypdot, asc, xinc, xp, yp, asctrm, v;
	double edot, cosf2, cose2, nudot, rdot, udot, idot, r_bds;
	double xsat0[3];
	double A0, AK, n0, delta_na;
	list<GPSEPH>::iterator gpsItr;
	GPSEPH* eph;
	intv = 2.0;
	switch (cprn[0]) {
	case 'G':
		GMEARTH = 3.986005E14;
		EARTH_ROTATE = 7.2921151467E-5;
		break;
	case 'E':
		GMEARTH = 3.986004418E14;
		EARTH_ROTATE = 7.2921151467E-5;
		break;
	case 'C':
		GMEARTH = 3.986004418E14;
		EARTH_ROTATE = 7.2921150E-5;
		break;
	case 'J':
		GMEARTH = 3.986005E14;
		EARTH_ROTATE = 7.2921150E-5;
		break;
	default:
		GMEARTH = GME;
		EARTH_ROTATE = E_ROTATE;
		break;
	}
	lfind = 0;
	*clk = 0;
	*dtmin = 4.0;
	isys = index_string(SYS, cprn[0]);
	double tmax = 0.0, ttag;
	for (gpsItr = this->gnssEph[isys].begin();
			gpsItr != this->gnssEph[isys].end(); ++gpsItr) {
		if (strstr((*gpsItr).cprn, cprn)) {
			dt = (wk - (*gpsItr).week) * 168.0 + (sow - (*gpsItr).toe) / 3600.0;
			if (*iode >= 0) {
				if(cprn[0] == 'C'){
					if ((*gpsItr).aode == *iode || (*gpsItr).iodc == *iode
							|| (*gpsItr).aodc == *iode) {
						*dtmin = fabs(dt);
						eph = &(*gpsItr);
						lfind = 1;
						break;
					}
				}else{
					if ((*gpsItr).aode == *iode ) {
						*dtmin = fabs(dt);
						eph = &(*gpsItr);
						lfind = 1;
						break;
					}
				}
			} else if (*iode == -1) {
				if (fabs(dt) <= *dtmin) {
					*dtmin = fabs(dt);
					eph = &(*gpsItr);
					lfind = 1;
				}
			} else if (*iode == -2) {
				ttag = (*gpsItr).mjd * 86400.0 + (*gpsItr).sod;
				// For the real-time process : always the newest
				if (ttag > tmax && fabs(dt) <= intv) {
					tmax = ttag;
					*dtmin = fabs(dt);
					lfind = 1;
					eph = &(*gpsItr);
				}
			}
		}
	}
	if (lfind == 0)
		return 0;
	*iode = eph->aode;
	dt = (wk - eph->week) * 604800.0 + sow - eph->toe;
	dtc = (wk * 7 + 44244 - eph->mjd) * 86400.0 + sow - eph->sod;

	*clk = eph->a0 + (eph->a1 + eph->a2 * dtc) * dtc;
	*tgd = eph->tgd;
	if(cprn[0] == 'C'){
		double f1 = BDS_B1;
		double f2 = BDS_B3;
		double a1 = f1 * f1 / (f1 * f1 - f2 * f2);
		double tgd =  eph->tgd; /**B1/B3  B2/B3*/
		*clk = *clk - a1 * tgd ;
	}
	if (strstr(mode, "ynn") != NULL)
		return 1;
	if (cprn[0] == 'C' && eph->signal_idx > 6) {
		f_prn = atoi(cprn + 1);
		if ((f_prn >= 38 && f_prn <= 40) || (f_prn >= 59 && f_prn <= 63)) {
			A0 = A_ref_igso + eph->delta_A;
		} else {
			A0 = A_ref_MEO + eph->delta_A;
		}
		a = AK = A0 + eph->A_DOT * dt;
		n0 = sqrt(GMEARTH / A0 / A0 / A0);
		delta_na = eph->dn + 0.5 * eph->delta_n_dot * dt;
		xn = n0 + delta_na;
		xm = eph->m0 + xn * dt;
	} else {
		a = eph->roota * eph->roota;
		xn = sqrt(GMEARTH / a / a / a);
		xn = xn + eph->dn;
		xm = eph->m0 + xn * dt;
	}
	ex = xm;
	e = eph->e;
	for (i = 0; i < 12; i++)
		ex = xm + e * sin(ex);
	v0 = 1.0 - e * cos(ex);
	vs = sqrt(1.0 - e * e) * sin(ex) / v0;
	vc = (cos(ex) - e) / v0;
	v = fabs(asin(vs));
	// decide vs
	if (vc >= 0) {
		if (vs < 0)
			v = 2.0 * PI - v;
	} else {
		if (vs <= 0)
			v = PI + v;
		else
			v = PI - v;
	}

	phi = v + eph->omega;
	ccc = cos(2.0 * phi);
	sss = sin(2.0 * phi);
	du = eph->cuc * ccc + eph->cus * sss;
	dr = eph->crc * ccc + eph->crs * sss;
	di = eph->cic * ccc + eph->cis * sss;
	r = a * (1.0 - e * cos(ex)) + dr;
	u = phi + du;

	xi = eph->i0 + eph->idot * dt + di;
	xx = r * cos(u);
	yy = r * sin(u);

	if (strstr(cprn, "C01") != NULL || strstr(cprn, "C02") != NULL
			|| strstr(cprn, "C03") != NULL || strstr(cprn, "C04") != NULL
			|| strstr(cprn, "C05") != NULL || strstr(cprn, "C59") != NULL) {
		xnode = eph->omega0 + eph->omegadot * dt;
	} else {
		xnode = eph->omega0 + (eph->omegadot - EARTH_ROTATE) * dt;
	}

	//IF THE TIME SYSTEM IS DIFFERENT,IT IS BETTER TO COMPUTE THE ORBIT IN THEIR OWN SYSTEM,
	//BECAUSE TOE WILL BE DIFFERENT IF YOU TRANSFORM IT TO OTHER SYSTEM

	xnode = xnode - EARTH_ROTATE * eph->toe;
	xsat[0] = xx * cos(xnode) - yy * cos(xi) * sin(xnode);
	xsat[1] = xx * sin(xnode) + yy * cos(xi) * cos(xnode);
	xsat[2] = yy * sin(xi);

	for (i = 0; i < 3; i++)
		svpos[i] = xsat[i];

	if (strstr(cprn, "C01") != NULL || strstr(cprn, "C02") != NULL
			|| strstr(cprn, "C03") != NULL || strstr(cprn, "C04") != NULL
			|| strstr(cprn, "C05") != NULL || strstr(cprn, "C59") != NULL) {
		// rx(-5.0*rad) for GEO
		pos[0] = xsat[0];
		pos[1] = cos(-5.0 * DEG2RAD) * xsat[1] + sin(-5.0 * DEG2RAD) * xsat[2];
		pos[2] = sin(5.0 * DEG2RAD) * xsat[1] + cos(-5.0 * DEG2RAD) * xsat[2];

		// rz(wearth*dt)
		xsat[0] = pos[0] * cos(EARTH_ROTATE * dt)
				+ pos[1] * sin(EARTH_ROTATE * dt);
		xsat[1] = -1.0 * pos[0] * sin(EARTH_ROTATE * dt)
				+ pos[1] * cos(EARTH_ROTATE * dt);
		xsat[2] = pos[2];
	}
	if (strstr(mode, "yyn") != NULL)
		return 1;
	asc = xnode;
	xinc = xi;
	xp = xx;
	yp = yy;

	edot = xn / (1 - e * cos(ex));
	cosf2 = cos(v / 2) * cos(v / 2);
	cose2 = cos(ex / 2) * cos(ex / 2);
	nudot = sqrt((1 + e) / (1 - e)) * cosf2 * edot / cose2;

	rdot = a * e * edot * sin(ex)
			+ 2 * (eph->crs * ccc - eph->crc * sss) * nudot;

	udot = (1 + 2.0 * eph->cus * ccc - 2.0 * eph->cuc * sss) * nudot;

	idot = eph->idot + 2.0 * (eph->cic * ccc - eph->cis * sss) * nudot;

	if (strstr(cprn, "C01") != NULL || strstr(cprn, "C02") != NULL
			|| strstr(cprn, "C03") != NULL || strstr(cprn, "C04") != NULL
			|| strstr(cprn, "C05") != NULL || strstr(cprn, "C59") != NULL ) {
		asctrm = eph->omegadot;
	} else {
		asctrm = eph->omegadot - EARTH_ROTATE;
	}

	xpdot = rdot * cos(u) - r * sin(u) * udot;
	ypdot = rdot * sin(u) + r * cos(u) * udot;

	xsat[3] = xpdot * cos(asc) - ypdot * cos(xinc) * sin(asc)
			+ yp * sin(asc) * sin(xinc) * idot - xp * sin(asc) * asctrm
			- yp * cos(xinc) * cos(asc) * asctrm;
	xsat[4] = xpdot * sin(asc) + ypdot * cos(xinc) * cos(asc)
			- yp * cos(asc) * sin(xinc) * idot + xp * cos(asc) * asctrm
			- yp * cos(xinc) * sin(asc) * asctrm;
	xsat[5] = ypdot * sin(xinc) + yp * cos(xinc) * idot;

	if (strstr(cprn, "C01") != NULL || strstr(cprn, "C02") != NULL
			|| strstr(cprn, "C03") != NULL || strstr(cprn, "C04") != NULL
			|| strstr(cprn, "C05") != NULL || strstr(cprn, "C59") != NULL) {
		// rx(-5.0*rad) for GEO
		geovel[0] = xsat[3];
		geovel[1] = cos(-5.0 * DEG2RAD) * xsat[4]
				+ sin(-5.0 * DEG2RAD) * xsat[5];
		geovel[2] = sin(5.0 * DEG2RAD) * xsat[4]
				+ cos(-5.0 * DEG2RAD) * xsat[5];

		// rz(wearth*dt)
		detavel[0] = geovel[0] * cos(EARTH_ROTATE * dt)
				+ geovel[1] * sin(EARTH_ROTATE * dt);
		detavel[1] = -1.0 * geovel[0] * sin(EARTH_ROTATE * dt)
				+ geovel[1] * cos(EARTH_ROTATE * dt);
		detavel[2] = geovel[2];

		// dot(R(w*t))*R(i)*XPOS
		geovel[0] = svpos[0];
		geovel[1] = cos(-5.0 * DEG2RAD) * svpos[1]
				+ sin(-5.0 * DEG2RAD) * svpos[2];
		geovel[2] = sin(5.0 * DEG2RAD) * svpos[1]
				+ cos(-5.0 * DEG2RAD) * svpos[2];

		// rz(wearth*dt)
		xsat[3] = EARTH_ROTATE
				* (-1.0 * geovel[0] * sin(EARTH_ROTATE * dt)
						+ geovel[1] * cos(EARTH_ROTATE * dt));
		xsat[4] = EARTH_ROTATE
				* (-1.0 * geovel[0] * cos(EARTH_ROTATE * dt)
						- geovel[1] * sin(EARTH_ROTATE * dt));
		xsat[5] = 0.0;

		// summary
		xsat[3] += detavel[0];
		xsat[4] += detavel[1];
		xsat[5] += detavel[2];
	}
	return 1;
}
int BroadcastEphUtils::m_gls2xyz(const char* mode, const char* cprn, int wk,
		double sow, double* xsat, double* clk, double* dtmin, int* iode) {
	double dt = 0.0, dt1 = 0.0, dt2 = 0.0, x[6], sod, intv;
	int i, nsign, lfind, mjd;
	list<GLSEPH>::iterator glsItr;
	GLSEPH* ephr;
	lfind = 0;
	*clk = 0.0;
	*dtmin = 12.0;
	wksow2mjd(wk, sow, &mjd, &sod);
	memset(x, 0, sizeof(double) * 6);
	double tmax = 0.0, ttag;
	intv = 2.0;
	for (glsItr = this->glsEph.begin(); glsItr != this->glsEph.end();
			++glsItr) {
		if (strstr((*glsItr).cprn, cprn)) {
			dt = (mjd - (*glsItr).mjd) * 24.0 + (sod - (*glsItr).sod) / 3600.0;
			if (*iode >= 0) {
				if (*iode == (*glsItr).aode) {
					*dtmin = fabs(dt);
					ephr = &(*glsItr);
					lfind = 1;
					break;
				}
			} else if (*iode == -1) {
				if (fabs(dt) <= *dtmin) {
					*dtmin = fabs(dt);
					ephr = &(*glsItr);
					lfind = 1;
				}
			} else if (*iode == -2) {
				ttag = (*glsItr).mjd * 86400.0 + (*glsItr).sod;
				if (ttag > tmax && fabs(dt) <= intv) {
					tmax = ttag;
					*dtmin = fabs(dt);
					ephr = &(*glsItr);
					lfind = 1;
				}
			}
		}
	}
	if (!lfind)
		return 0;
	*iode = ephr->aode;
	dt = (mjd - ephr->mjd) * 86400.0 + sod - ephr->sod;
	*clk = ephr->tau;
	if (strstr(mode, "ynn") != NULL)
		return 1;

	dt1 = dt - (int) (dt / 60.0) * 60.0;
	dt2 = dt1;
	nsign = 1;
	if (dt < 0.0)
		nsign = -1;

	glsinit(x, *ephr);

	if (dt1 != 0.0)
		glsrkf4(dt1, x, *ephr);

	while (ABS(dt - dt1) >= 60.0) {
		if (dt - dt1 != 0.0) {
			glsrkf4(nsign * 60.0, x, *ephr);
			dt1 = dt1 + nsign * 60.0;
		}
	}
	pz902wgs84(mjd, sod, x, xsat, "MCC");
	for (i = 0; i < 3; i++)
		xsat[i] = xsat[i] * 1e3;
	for (i = 3; i < 6; i++)
		xsat[i] = x[i] * 1e3;
	return 1;
}

void BroadcastEphUtils::glsinit(double* x, GLSEPH& eph) {
	int i;
	for (i = 0; i < 3; i++)
		x[i] = eph.pos[i];
	for (i = 3; i < 6; i++)
		x[i] = eph.vel[i - 3];
}
void BroadcastEphUtils::glsrkf4(double h, double* x, GLSEPH& eph) {
	int i, j;
	double coef[4], xj[6], acc[6], funct[6];
	coef[0] = h / 2;
	coef[1] = h / 2;
	coef[2] = h;
	coef[3] = h;
	for (i = 0; i < 6; i++)
		xj[i] = x[i];
	glsfright(xj, acc, eph);
	for (i = 0; i < 6; i++)
		funct[i] = xj[i];
	for (j = 0; j < 3; j++) {
		for (i = 0; i < 6; i++) {
			xj[i] = x[i] + coef[j] * acc[i];
			funct[i] = funct[i] + coef[j + 1] * acc[i] / 3.0;
		}
		glsfright(xj, acc, eph);
	}
	for (i = 0; i < 6; i++)
		x[i] = funct[i] + h * acc[i] / 6.0;
}
void BroadcastEphUtils::glsfright(double* x, double* acc, GLSEPH& eph) {
	const double c20 = -1.08262575e-3;
	const double GMEARTH = 3.986004418e5;
	const double EROT = 7.292115e-5;
	const double ER = 6.378136e3;
	int i;
	double factor1, factor2, r, vec1[3], vec2[3];
	r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	factor1 = -GMEARTH / pow(r, 3);
	factor2 = 1.5 * GMEARTH * c20 * ER * ER / pow(r, 5);
	vec1[0] = EROT * EROT * x[0] + 2.0 * EROT * x[4];
	vec1[1] = EROT * EROT * x[1] - 2.0 * EROT * x[3];
	vec1[2] = 0.0;
	vec2[0] = x[0] * (1.0 - 5.0 * pow((x[2] / r), 2));
	vec2[1] = x[1] * (1.0 - 5.0 * pow((x[2] / r), 2));
	vec2[2] = x[2] * (3.0 - 5.0 * pow((x[2] / r), 2));
	for (i = 0; i < 3; i++) {
		acc[i] = x[i + 3];
		acc[i + 3] = eph.acc[i] + vec1[i] + factor2 * vec2[i] + factor1 * x[i];
	}
}
void BroadcastEphUtils::pz902wgs84(int mjd, double sod, double* pos,
		double *xsat, const char* trans) {
	double factor, biase[3], mat[3][3], vec[3];
	int i, j;
	if (timdif(mjd, sod, 54362, 86367.0) >= 0.0) {
		factor = 1.0;
//		biase[0] = -0.36e-3;
//		biase[1] = 0.08e-3;
//		biase[2] = 0.18e-3;

		biase[0] = 0.008e-3;
		biase[1] = 0.001e-3;
		biase[2] = 0.001e-3;
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++) {
				if (i == j) {
					mat[i][j] = 1.0;
				} else {
					mat[i][j] = 0.0;
				}
			}
	} else {
		if (strstr(trans, "MCC") != NULL) {
			factor = 1.0 + 22.0e-9;
			biase[0] = -0.47e-3;
			biase[1] = -0.51e-3;
			biase[2] = -1.56e-3;
			mat[0][0] = 1.0;
			mat[0][1] = -1.728e-6;
			mat[0][2] = -1.7e-8;
			mat[1][0] = 1.728e-6;
			mat[1][1] = 1.0;
			mat[1][2] = 7.6e-8;
			mat[2][0] = 1.7e-8;
			mat[2][1] = -7.6e-8;
			mat[2][2] = 1.0;
		} else if (strstr(trans, "RUS") != NULL) {
			factor = 1.0 - 1.2e-7;
			biase[0] = -1.1e-3;
			biase[1] = -0.3e-3;
			biase[2] = -0.9e-3;
			mat[0][0] = 1.0;
			mat[0][1] = -8.2e-7;
			mat[0][2] = 0.0;
			mat[1][0] = 8.2e-7;
			mat[1][1] = 1.0;
			mat[1][2] = 0.0;
			mat[2][0] = 0.0;
			mat[2][1] = 0.0;
			mat[2][2] = 1.0;
		} else if (strstr(trans, "LMU") != NULL) {
			factor = 1.0;
			biase[0] = 0.0;
			biase[1] = 0.0;
			biase[2] = 0.0;
			mat[0][0] = 1.0;
			mat[0][1] = -1.6e-6;
			mat[0][2] = 0.0;
			mat[1][0] = 1.6e-6;
			mat[1][1] = 1.0;
			mat[1][2] = 0.0;
			mat[2][0] = 0.0;
			mat[2][1] = 0.0;
			mat[2][2] = 1.0;
		} else if (strstr(trans, "MIT") != NULL) {
			factor = 1.0;
			biase[0] = 0.0;
			biase[1] = 2.5e-3;
			biase[2] = 0.0;
			mat[0][0] = 1.0;
			mat[0][1] = -1.9e-6;
			mat[0][2] = 0.0;
			mat[1][0] = 1.9e-6;
			mat[1][1] = 1.0;
			mat[1][2] = 0.0;
			mat[2][0] = 0.0;
			mat[2][1] = 0.0;
			mat[2][2] = 1.0;
		} else {
			printf("***ERROR(pz902wgs84): unknown transformation type %s!\n",
					trans);
			exit(1);
		}
	}
	matmpy((double*) mat, pos, vec, 3, 3, 1);
	for (i = 0; i < 3; i++)
		xsat[i] = biase[i] + factor * vec[i];
}

void RnxNavFile::m_readRnxnav(char csys) {
	// read rinex navigation file
	int i, k, isys, len, iy, im, id, ih, imin, mjd, wk_in;
	double dsec, dt1, dt0, d_aode, sod;
	char buf[2048], cline[1024], *ptr;
	list<GPSEPH>::iterator gpsItr;
	list<GLSEPH>::iterator glsItr;
	GPSEPH eph;
	GLSEPH ephr;
	ifstream in;
	string line;
	streampos pos;
	bool already;
	double data1, data2, data3;	//for BDS B2b
	in.open(curFile, ios::in | ios::binary);
	if (!in.is_open()) {
		cout << "***ERROR(m_readRnxnav):can't open file " << curFile << endl;
		exit(1);
	}
	while (getline(in, line)) {
		if (strstr(line.c_str(), "END OF HEADER"))
			break;
		if (strstr(line.c_str(), "RINEX VERSION / TYPE") != NULL) {
			sscanf(line.c_str(), "%lf", &ver);
			if (ver > 3.0) {
				if (csys == 'M' && line[40] != 'M') {
					printf(
							"###WARNING(read_rnxnav):there only broadcast for %c!\n",
							line[40]);
				}
			}
		} else if (strstr(line.c_str(), "PGM / RUN BY / DATE") != NULL)
			continue;
		else if (strstr(line.c_str(), "COMMENT") != NULL)
			continue;
		else if (strstr(line.c_str(), "ION ALPHA") != NULL) {

		} else if (strstr(line.c_str(), "ION BETA") != NULL) {

		}
		//rinex 3.00 3.01 3.02
		else if (strstr(line.c_str(), "IONOSPHERIC CORR") != NULL) {
			if (!strncmp(line.c_str(), "GPS ", 4)
					|| !strncmp(line.c_str(), "GPSA", 4)) {
				i = index_string(SYS, 'G');
				k = 0;
			} else if (!strncmp(line.c_str(), "GPSB", 4)) {
				i = index_string(SYS, 'G');
				k = 1;
			} else if (!strncmp(line.c_str(), "GAL ", 4)) {
				i = index_string(SYS, 'E');
				k = 0;
			} else if (!strncmp(line.c_str(), "BDS ", 4)
					|| !strncmp(line.c_str(), "BDSA", 4)) {
				i = index_string(SYS, 'C');
				k = 0;
			} else if (!strncmp(line.c_str(), "BDSB", 4)) {
				i = index_string(SYS, 'C');
				k = 1;
			} else if (!strncmp(line.c_str(), "QZS ", 4)
					|| !strncmp(line.c_str(), "QZSA", 4)) {
				i = index_string(SYS, 'J');
				k = 0;
			} else if (!strncmp(line.c_str(), "QZSB", 4)) {
				i = index_string(SYS, 'J');
				k = 1;
			} else
				printf("###WARNING(read_rnxnav):unknown ionospheric corr!\n");
			strncpy(ionc[i][k], line.c_str(), 4);
			strcpy(buf, line.c_str());
			ptr = buf;
			while (*ptr != '\0') {
				if (*ptr == 'D')
					*ptr = 'e';
				ptr++;
			}
			sscanf(buf, "%*s%lf%lf%lf%lf", ion[i][k], ion[i][k] + 1,
					ion[i][k] + 2, ion[i][k] + 3);
		} else if (strstr(line.c_str(), "DELTA-UTC: A0,A1,T,W") != NULL) {
			i = index_string(SYS, csys);
			if (i == -1) {
				printf(
						"$$$MESSAGE(read_rnxnav):DELTA-UTC: A0,A1,T,W is only valid for single system!\n");
				continue;
			}
			strcpy(buf, line.c_str());
			ptr = buf;
			while (*ptr != '\0') {
				if (*ptr == 'D')
					*ptr = 'e';
				ptr++;
			}
			sscanf(buf, "%lf%lf%lf%lf", tim[i][0], tim[i][0] + 1, tim[i][0] + 2,
					tim[i][0] + 3);
		} else if (strstr(line.c_str(), "TIME SYSTEM CORR") != NULL) {
			if (!strncmp(line.c_str(), "GPUT", 4)) {
				i = index_string(SYS, 'G');
				k = 0;
			} else if (!strncmp(line.c_str(), "GLUT", 4)) {
				i = index_string(SYS, 'R');
				k = 0;
			} else if (!strncmp(line.c_str(), "GAUT", 4)) {
				i = index_string(SYS, 'E');
				k = 0;
			} else if (!strncmp(line.c_str(), "BDUT", 4)) {
				i = index_string(SYS, 'C');
				k = 0;
			} else if (!strncmp(line.c_str(), "SBUT", 4)) {
				i = index_string(SYS, 'S');
				k = 0;
			} else if (!strncmp(line.c_str(), "QZUT", 4)) {
				i = index_string(SYS, 'J');
				k = 0;
			} else if (!strncmp(line.c_str(), "GAGP", 4)) {
				i = index_string(SYS, 'G');
				k = 1;
			} else if (!strncmp(line.c_str(), "GLGP", 4)) {
				i = index_string(SYS, 'R');
				k = 1;
			} else if (!strncmp(line.c_str(), "QZGP", 4)) {
				i = index_string(SYS, 'J');
				k = 1;
			} else {
				printf(
						"$$$MESSAGE(read_rnxnav):unknown unknown TIME SYSTEM CORR!\n");
				continue;
			}
			strncpy(timc[i][k], line.c_str(), 4);
			strcpy(buf, line.c_str());
			ptr = buf;
			while (*ptr != '\0') {
				if (*ptr == 'D')
					*ptr = 'e';
				ptr++;
			}
			sscanf(buf, "%*5c%lf%lf%lf%lf", tim[i][k], tim[i][k] + 1,
					tim[i][k] + 2, tim[i][k] + 3);
		} else if (strstr(line.c_str(), "LEAP SECONDS") != NULL)
			sscanf(line.c_str(), "%d", &leap);
	}
	while (!in.eof()) {
		// acquire the ephemeris
		pos = in.tellg();
		getline(in, line);
		if (!len_trim(line.c_str()))
			continue;
		in.seekg(pos, ios::beg);
		if (ver < 3.0) {
			buf[0] = csys;
		} else
			buf[0] = line[0];
		switch (buf[0]) {
		case 'G':
		case 'E':
		case 'C':
		case 'J':
			isys = index_string(SYS, buf[0]);
			len = 0;
			memset(buf, 0, sizeof(char) * 2048);
			for (i = 0; i < 8; i++) {
				getline(in, line);
				strcpy(cline, line.c_str());
				filleph(cline, ver);
				/*if (i != 8) {
				 filleph(cline, ver);
				 }*/
				if (ver < 3) {
					if (buf[0] != csys)
						buf[0] = csys;
					strncpy(buf + 1 + len, cline, strlen(cline));
				} else
					strncpy(buf + len, cline, strlen(cline));
				len += strlen(cline);
			}
			ptr = buf;
			while (*ptr != '\0') {
				if (*ptr == '\n' || *ptr == '\r')
					*ptr = ' ';

				if (*ptr == 'D')
					*ptr = 'e';
				ptr++;
			}
			if (ver < 3.0) {
				buf[0] = csys;
				if (buf[1] == ' ')
					buf[1] = '0';
			}
			if (buf[0] == 'C') {	//BDS B2b
				sscanf(buf,
						"%s%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
						eph.cprn, &iy, &im, &id, &ih, &imin, &dsec, &eph.a0,
						&eph.a1, &eph.a2, &d_aode, &eph.crs, &eph.dn, &eph.m0,
						&eph.cuc, &eph.e, &eph.cus, &eph.roota, &eph.toe,
						&eph.cic, &eph.omega0, &eph.cis, &eph.i0, &eph.crc,
						&eph.omega, &eph.omegadot, &eph.idot, &eph.resvd0,
						&eph.week, &eph.resvd1, &eph.accu, &eph.hlth, &eph.tgd,
						&eph.aodc, &data1, &eph.iodc, &data3, &eph.signal_idx);
				eph.signal_idx = int(eph.signal_idx) == 0 ? 1 : int(eph.signal_idx);
				if (eph.signal_idx > 6 && eph.signal_idx <= 13) {
					eph.delta_A = eph.roota;
					eph.A_DOT = eph.resvd0;
					eph.delta_n_dot = eph.resvd1;
				}
				if (eph.signal_idx == 1 || eph.signal_idx == 0) { //B1I
					eph.tgd_BDS[0] = eph.tgd;
					eph.isc_BDS[0] = 0.0;
				} else if (eph.signal_idx == 2) { //B2I
					eph.tgd_BDS[1] = eph.aodc;
					eph.isc_BDS[1] = 0.0;
				} else if (eph.signal_idx == 3) { //B3I
					eph.tgd_BDS[2] = 0.0;
					eph.isc_BDS[2] = 0.0;
				} else if (eph.signal_idx == 4) { //B1Q
					eph.tgd_BDS[3] = eph.tgd;
					eph.isc_BDS[3] = 0.0;
				} else if (eph.signal_idx == 5) { //B1Q
					eph.tgd_BDS[4] = eph.aodc;
					eph.isc_BDS[4] = 0.0;
				} else if (eph.signal_idx == 6) { //B3Q
					eph.tgd_BDS[5] = 0.0;
					eph.isc_BDS[5] = 0.0;
				} else if (eph.signal_idx == 7) { //B1C
					eph.tgd_BDS[6] = eph.tgd;
					eph.isc_BDS[6] = data3;
				} else if (eph.signal_idx == 8) { //B2a
					eph.tgd_BDS[7] = eph.tgd;
					eph.isc_BDS[7] = data3;
				} else if (eph.signal_idx == 9) { //B2bI
					eph.tgd_BDS[8] = eph.tgd;
					eph.isc_BDS[8] = 0.0;
				} else if (eph.signal_idx == 10) { //B2bQ
					eph.tgd_BDS[9] = eph.tgd;
					eph.isc_BDS[9] = 0.0;
				} else if (eph.signal_idx == 11) { //B1A
					eph.tgd_BDS[10] = eph.tgd;
					eph.isc_BDS[10] = data3;
				} else if (eph.signal_idx == 12) { //B3A
					eph.tgd_BDS[11] = 0.0;
					eph.isc_BDS[11] = data3;
				} else if (eph.signal_idx == 13) { //B3AE
					eph.tgd_BDS[12] = eph.aodc;
					eph.isc_BDS[12] = data3;
				}
//				if(eph.signal_idx != 1 && eph.signal_idx != 0) /**skip if that is not 18 parameters **/
				if(eph.signal_idx != 7) /**skip if that is not 18 parameters **/
					continue;

			} else {
				sscanf(buf,
						"%s%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
						eph.cprn, &iy, &im, &id, &ih, &imin, &dsec, &eph.a0,
						&eph.a1, &eph.a2, &d_aode, &eph.crs, &eph.dn, &eph.m0,
						&eph.cuc, &eph.e, &eph.cus, &eph.roota, &eph.toe,
						&eph.cic, &eph.omega0, &eph.cis, &eph.i0, &eph.crc,
						&eph.omega, &eph.omegadot, &eph.idot, &eph.resvd0,
						&eph.week, &eph.resvd1, &eph.accu, &eph.hlth, &eph.tgd,
						&eph.aodc, &data1, &eph.iodc, &data3, &eph.signal_idx);
			}
			/////////////////////////// for BDS III test,must be removed //////////////////////////////////
			if (atoi(eph.cprn + 1) > 14) {
				eph.hlth = 0.0;
			}
			///////////////////////////////////////////////////////////////////////////////////////////////
			if (eph.hlth > 0 || isys == -1)
				continue;
			yr2year(iy);
			if (iy < 2000)
				continue;
			// time of clock
			eph.mjd = md_julday(iy, im, id); // TIME IN BDS TIME AND SHOULD NOT CHANGE IT INTO GPST BECAUSE THERE ARE OTHER TIME TAG IN THE EPHEMERIS
			eph.sod = ih * 3600.0 + imin * 60.0 + dsec;

			// adapter to bds
			if (eph.cprn[0] == 'C') {
				//the week in broadcast file generated by WHU is GPS week,
				//but that in IGS meraged file is BDS week
				mjd2wksow(eph.mjd, eph.sod, &k, &dt1);
				if (eph.week != 0.0) {
					wksow2mjd(eph.week, eph.toe, &mjd, &sod);
					mjd2wksow(mjd, sod, &wk_in, &eph.toe);
					eph.week = wk_in;
					if (k != (int) eph.week) {
						eph.week = eph.week + 1356.0;
					}
				} else {
					eph.week = 0;
					for (int i = -1; i < 2; i++) {
						wksow2mjd(k + i, eph.toe, &mjd, &sod);
						if (fabs(timdif(mjd, sod, eph.mjd, eph.sod))
								< 86400.0 * 7 / 2) {
							eph.week = k + i;
							break;
						}
					}
					if (eph.week == 0) {
						memset(eph.cprn, 0, sizeof(eph.cprn));
					} else {
						wksow2mjd(eph.week, eph.toe, &mjd, &sod);
						mjd2wksow(mjd, sod, &wk_in, &eph.toe);
						eph.week = wk_in;
					}
				}
				eph.tgd1 = eph.aodc;
			}
			//check time
			dt0 = 0.0;
			dt1 = 0.0;
			if (mjd0 != 0) {
				dt0 = eph.mjd + eph.sod / 86400.0 - mjd0;
			}
			if (mjd1 != 0) {
				dt1 = eph.mjd + eph.sod / 86400.0 - mjd1;
			}
			if (dt0 < -1.0 / 24.0 || dt1 > 1.0 / 24.0)
				continue;
			already = false;
			for (gpsItr = this->gnssEph[isys].begin();
					gpsItr != this->gnssEph[isys].end(); gpsItr++) {
				if (strstr((*gpsItr).cprn, eph.cprn) && (*gpsItr).mjd == eph.mjd
						&& (*gpsItr).sod == eph.sod) {
					already = true;
					break;
				}
			}
			if (!already) {
				neph[isys] = neph[isys] + 1;
				eph.aode = genAode(buf[0], eph.mjd, eph.sod, eph.toe,
						static_cast<int>(d_aode), &eph);
				this->gnssEph[isys].push_back(eph);
			}
			break;
		case 'R':
			isys = index_string(SYS, buf[0]);
			len = 0;
			memset(buf, 0, sizeof(char) * 2048);
			for (i = 0; i < 4; i++) {
				getline(in, line);
				strcpy(cline, line.c_str());
				if (ver < 3) {
					if (buf[0] != csys)
						buf[0] = csys;
					strncpy(buf + 1 + len, cline, strlen(cline));
				} else
					strncpy(buf + len, cline, strlen(cline));
				len += strlen(cline);
			}
			ptr = buf;
			while (*ptr != '\0') {
				if (*ptr == '\n' || *ptr == '\r') {
					*ptr = ' ';
				}
				if (*ptr == 'D') {
					*ptr = 'e';
				}
				ptr++;
			}
			if (ver < 3.0) {
				buf[0] = csys;
				if (buf[1] == ' ')
					buf[1] = '0';
			}
			sscanf(buf,
					"%s%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
					ephr.cprn, &iy, &im, &id, &ih, &imin, &dsec, &ephr.tau,
					&ephr.gamma, &ephr.tk, &ephr.pos[0], &ephr.vel[0],
					&ephr.acc[0], &ephr.health, &ephr.pos[1], &ephr.vel[1],
					&ephr.acc[1], &ephr.frenum, &ephr.pos[2], &ephr.vel[2],
					&ephr.acc[2], &ephr.age);

			if (ephr.health > 0.0)
				continue;
			yr2year(iy);
			if (iy < 2000)
				continue;
			ephr.mjd = md_julday(iy, im, id);
			ephr.sod = ih * 3600.0 + imin * 60.0 + dsec;
			brdtime(ephr.cprn, &ephr.mjd, &ephr.sod); // CHANGE IT INTO GPST TIME BECAUSE THE THE BASE TIME CAN CHANGE AND SHOULD NOT
													  // UPDATE IT IN A WHILE LOOP
			//check time
			dt0 = 0.0;
			dt1 = 0.0;
			if (mjd0 != 0) {
				dt0 = ephr.mjd + ephr.sod / 86400.0 - mjd0;
			}
			if (mjd1 != 0) {
				dt1 = ephr.mjd + ephr.sod / 86400.0 - mjd1;
			}
			if (dt0 < -3.0 / 24.0 || dt1 > 3.0 / 24.0)
				continue;
			already = false;
			for (glsItr = this->glsEph.begin(); glsItr != this->glsEph.end();
					glsItr++) {
				if (strstr((*glsItr).cprn, ephr.cprn)
						&& (*glsItr).mjd == ephr.mjd
						&& (*glsItr).sod == ephr.sod) {
					already = true;
					break;
				}
			}
			if (!already) {
				neph[isys] = neph[isys] + 1;
				ephr.aode = genAode(buf[0], ephr.mjd, ephr.sod, 0.0, ephr.aode,
				NULL);
				this->glsEph.push_back(ephr);
			}
			break;
		default:
			getline(in, line);
			break;
		}
	}
	in.close();
	cout << "The epoches of GPS is:" << this->gnssEph[0].size() << endl;
	cout << "The epoches of BDS is:" << this->gnssEph[3].size() << endl;
}

int RnxNavFile::v_readOrbit(const char* cprn, int mjd, double sod,
		double* xsat) {
	// using rinex navigation file to acquire orbit
	int wk, iode = -1;
	double sow, tgd, dtmin, sclk;
	mjd2wksow(mjd, sod, &wk, &sow);
	if (cprn[0] != 'R')
		if (cprn[0] == 'C')
			return this->m_brd2xyz("yyy", cprn, wk, sow - 14.0, xsat, &sclk,
					&dtmin, &tgd, &iode);
		else
			return this->m_brd2xyz("yyy", cprn, wk, sow, xsat, &sclk, &dtmin,
					&tgd, &iode);
	else
		return this->m_gls2xyz("yyy", cprn, wk, sow, xsat, &sclk, &dtmin, &iode);
}

int RnxNavFile::v_readClk(const char* cprn, int mjd, double sod, double* sclk) {
	// using rinex navigation file to acquire satellite clock
	int wk, iode = -1;
	double sow, tgd, dtmin;
	mjd2wksow(mjd, sod, &wk, &sow);
	*sclk = 0.0;
	if (cprn[0] != 'R')
		if (cprn[0] == 'C')
			return this->m_brd2xyz("ynn", cprn, wk, sow - 14.0, NULL, sclk,
					&dtmin, &tgd, &iode);
		else
			return this->m_brd2xyz("ynn", cprn, wk, sow, NULL, sclk, &dtmin,
					&tgd, &iode);
	else
		return this->m_gls2xyz("ynn", cprn, wk, sow, NULL, sclk, &dtmin, &iode);
}

void RnxNavFile::v_openRnxEph(string ephcmd) {
	int curmjd;
	this->isOpen = 1;
	int iyear, imon, iday, ih, imin;
	double dsec, cursod;
	curmjd = atoi(ephcmd.substr(0, index_string(ephcmd.c_str(), ':')).c_str());
	cursod = atof(ephcmd.substr(index_string(ephcmd.c_str(), ':') + 1).c_str());
	mjd0 = curmjd;
	sod0 = 0.0;
	mjd1 = curmjd + 1;
	sod1 = 0.0;
	mjd2date(curmjd, cursod, &iyear, &imon, &iday, &ih, &imin, &dsec);
	memset(curFile, 0, sizeof(char) * 1024);
	for (int i = 0; i < this->navfilename.length(); i++)
		this->curFile[i] = this->navfilename[i];

	/**** using rinex3.x format to read nav ****/
	this->m_readRnxnav('M');
}

int RnxNavFile::v_readOrbit(const char* cprn, int mjd, double sod, double* xsat,
		int* iode) {
	// using rinex navigation file to acquire orbit
	// using rinex navigation file to acquire orbit
	int wk;
	double sow, tgd, dtmin, sclk;
	mjd2wksow(mjd, sod, &wk, &sow);
	if (cprn[0] != 'R')
		if (cprn[0] == 'C')
			return this->m_brd2xyz("yyy", cprn, wk, sow - 14.0, xsat, &sclk,
					&dtmin, &tgd, iode);
		else
			return this->m_brd2xyz("yyy", cprn, wk, sow, xsat, &sclk, &dtmin,
					&tgd, iode);
	else
		return this->m_gls2xyz("yyy", cprn, wk, sow, xsat, &sclk, &dtmin, iode);
}

int RnxNavFile::v_readClk(const char* cprn, int mjd, double sod, double* sclk,
		int* iode) {
	int wk;
	double sow, tgd, dtmin;
	mjd2wksow(mjd, sod, &wk, &sow);
	*sclk = 0.0;
	if (cprn[0] != 'R')
		if (cprn[0] == 'C')
			return this->m_brd2xyz("ynn", cprn, wk, sow - 14.0, NULL, sclk,
					&dtmin, &tgd, iode);
		else
			return this->m_brd2xyz("ynn", cprn, wk, sow, NULL, sclk, &dtmin,
					&tgd, iode);
	else
		return this->m_gls2xyz("ynn", cprn, wk, sow, NULL, sclk, &dtmin, iode);
}
