/*
 * Broadcast.h
 *
 *  Created on: 2020年5月11日
 *      Author: dongdong
 */

#ifndef INCLUDE_BROADCAST_H_
#define INCLUDE_BROADCAST_H_
#include <vector>
#include <list>
#include <string>
#include "Const.h"

using namespace std;

class ORBHD {
public:
	ORBHD() {
		nprn = nequ = rmjd = mjd0 = mjd1 = 0;
		sod0 = sod1 = rsod = dintv = 0;
		memset(cprn, 0, sizeof(cprn));
	}
	int nprn;
	char cprn[MAXSAT][LEN_PRN];
	int mjd0;
	int mjd1;
	int rmjd;
	int nequ;
	double sod0;
	double sod1;
	double rsod;
	double dintv;
};

class GPSEPH {
public:
	GPSEPH() {
		mjd = 0;
		sod = 0.0;
		iodc=-1;
		signal_idx=0;
		delta_A=0.0;
		A_DOT=0.0;
		delta_n_dot=0.0;
		memset(cprn, 0, sizeof(cprn));
		time(&gent);
	}
	char cprn[LEN_PRN];
	int mjd;
	double sod;
	// a[0]: SV clock offset
	// a[1]: SV clock drift
	// a[2]: SV clock drift rate
	double a0, a1, a2;
	// aode: age of ephemeris upload
	// crs, crc: Ortital radius correction
	// dn: Mean motion difference
	// m0: Mean anomaly at reference epoch
	// e: Eccentricity
	// cuc, cus: Latitude argument correction
	// roota: Square root of semi-major axis
	unsigned int aode;
	double crs, dn;
	double m0, cuc, e;
	double cus, roota;
	// toe, week: Ephemerides reference epoch in seconds with the week
	// cis, cic: Inclination correction
	// omega0: Longtitude of ascending node at the begining of the week
	// i0: Inclination at reference epoch
	// omega: Argument of perigee
	// omegadot: Rate of node's right ascension
	double toe, cic, omega0;
	double cis, i0, crc;
	double omega, omegadot;
	// idot: Rate of inclination angle
	// sesvd0:
	// resvd1:
	// accu: SV accuracy
	// hlth: SV health
	// tgd: Time group delay
	// aodc: Age of clock parameter upload
	double idot, resvd0, week, resvd1;
	double accu, hlth, tgd, tgd1, aodc;
	time_t gent;
	//for BDS-B2b
	double signal_idx;
	double iodc;
	double tgd_BDS[13], isc_BDS[13];
	double delta_A, A_DOT, delta_n_dot;
};

class GLSEPH {
public:
	GLSEPH() {
		mjd = 0;
		sod = 0.0;
		memset(cprn, 0, sizeof(cprn));
		time(&gent);
	}
	char cprn[LEN_PRN];
	int mjd, aode;
	double sod;
	// tau: SV clock bias
	// gama: SV relative frequency bias
	// tk: message frame time (tk+nd*86400) in seconds of the UTC week
	// pos: coordinate at ephemerides reference epoch in PZ-90
	// vel: velocity at ephemerides reference epoch in PZ-90
	// acc: acceleration at ephemerides reference epoch in PZ-90
	// health: SV health
	// frenum: frequency number
	// age: age of operation information
	double tau;
	double gamma;
	double tk;
	double pos[3];
	double vel[3];
	double acc[3];
	double health;
	double frenum;
	double age;
	time_t gent;
};

class ORBCLK {
public:
	ORBCLK() {
		int isat;
		wk = 0;
		sow = 0.0;
		memset(cprn, 0, sizeof(char) * MAXSAT * LEN_PRN);
		memset(c, 0, sizeof(double) * MAXSAT * 3);
		memset(dx, 0, sizeof(double) * MAXSAT * 6);
		for (isat = 0; isat < MAXSAT; isat++)
			iod[isat] = -1;
	}
	int wk;
	double sow;
	char cprn[MAXSAT][LEN_PRN];
	int iod[MAXSAT];
	double c[MAXSAT][3];
	double dx[MAXSAT][6];
};

class BroadcastEphUtils {
public:
	BroadcastEphUtils() {
		memset(neph, 0, sizeof(neph));
	}
	virtual ~BroadcastEphUtils() {
	}
	virtual int m_brd2xyz(const char* mode, const char* cprn, int wk,
			double sow, double* xsat, double* clk, double* dtmin, double* tgd,
			int* iode);
	virtual int m_gls2xyz(const char* mode, const char* cprn, int wk,
			double sow, double* xsat, double* clk, double* dtmin, int* iode);
	void glsinit(double* x, GLSEPH& eph);
	void glsrkf4(double h, double* x, GLSEPH& eph);
	void glsfright(double* x, double* acc, GLSEPH& eph);
	void pz902wgs84(int mjd, double sod, double* pos, double *xsat,
			const char* trans);
	int neph[MAXSYS];
	list<GPSEPH> gnssEph[MAXSYS];
	list<GLSEPH> glsEph;
};

class RnxNavFile: public BroadcastEphUtils {
public:
	virtual void v_openRnxEph(string ephcmd);
	void m_readRnxnav(char csys);
	virtual int v_readOrbit(const char* cprn, int mjd, double sod,
			double* xsat);
	virtual int v_readClk(const char* cprn, int mjd, double sod, double* sclk);
	virtual int v_readOrbit(const char* cprn, int mjd, double sod, double* xsat,
			int* iode);
	virtual int v_readClk(const char* cprn, int mjd, double sod, double* sclk,
			int* iode);
	char curFile[1024];
	int mjd0, mjd1;
	double sod0, sod1;
	int isOpen;
	char ionc[MAXSYS][2][LEN_EPHION];
	double ion[MAXSYS][2][4];
	string navfilename;
private:
	// header part
	double ver;
	// corrections to transform the system to UTC or other time systems
	char timc[MAXSYS][2][LEN_EPHION];
	double tim[MAXSYS][2][4];
	// number of leap second since 6-Jan-980
	int leap;
};
#endif /* INCLUDE_BROADCAST_H_ */
