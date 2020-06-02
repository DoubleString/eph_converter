/*
 * B2bDecoder.h
 *
 *  Created on: 2020年5月12日
 *      Author: dongdong
 */

#ifndef INCLUDE_B2BDECODER_H_
#define INCLUDE_B2BDECODER_H_
#include <vector>
#include "Const.h"
#include "Com.h"
#include "Broadcast.h"

using namespace std;

class sat_orb {
public:
	char cprn[LEN_PRN];
	int iode;
	double orb[6];
};

class sat_clk {
public:
	char cprn[LEN_PRN];
	int iode;
	double clk[3];
};
//orbit for one epoch
class orbit {
public:
	int mjd;
	double sod;
	map<string,sat_orb> ORB;
};
//clk for one epoch
class clk {
public:
	int mjd;
	double sod;
	map<string,sat_clk> CLK;
};

class B2b {
public:
	B2b();
	void m_process();
	int mjd0,mjd1,lbdtime;
	double sod0,sod1;
	RnxNavFile nav;
	char File[1024];
protected:
	void m_geneph_brd();
	void m_genclock_brd();
	void start_decode();
	void m_loadorbclk();
	void m_geneph();
	void m_genclock();
	void m_wtclkheader(FILE* fp);
	void m_wtephheader(FILE* fp);
	vector<orbit> ORB_ALL;
	vector<clk> CLK_ALL;
	vector<string> cprn_list;
};

#endif /* INCLUDE_B2BDECODER_H_ */
