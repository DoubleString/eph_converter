/*
 * B2bDecoder.cpp
 *
 *  Created on: 2020年5月12日
 *      Author: juntao
 */
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <fstream>
#include <vector>
#include "stdio.h"
#include "../include/B2bDecoder.h"

using namespace std;
B2b::B2b(){
	lbdtime = false;
	mjd0 = mjd1 = 0;
	sod0 = sod1 = 0.0;
	/******* if want to use the whole prn list,use this code to update cprn_list,other wise will use  ******/
	cprn_list.clear();
	for(int i = 0 ; i < 32;i++){
		char cprn[4] = {0};
		sprintf(cprn,"G%02d",i + 1);
		cprn_list.push_back(cprn);
	}
	for(int i = 0 ; i < 59;i++){
		char cprn[4] = {0};
		sprintf(cprn,"C%02d",i + 1);
		cprn_list.push_back(cprn);
	}
	/***********************************************************************************/
}
void B2b::m_loadorbclk(){
	string line;
	ifstream in;
	int iyear,imon,iday,ih,imin,iv,nsat;
	double dsec;
	in.open(File, ios::in | ios::binary);
	if (!in.is_open()) {
		cout << "***ERROR(m_loadorbclk):can't open file " << File << endl;
		exit(1);
	}
	cout << "start reading file..." << endl;
	while (!in.eof()) {
		getline(in, line);
		if (!len_trim(line.c_str()))
			continue;
		if (strstr(line.c_str(), "> ORBIT") != nullptr) {
			orbit orb_epoch;
			sscanf(line.c_str(),"%*s%*s%d%d%d%d%d%lf%d%d",&iyear,&imon,&iday,&ih,&imin,&dsec,&iv,&nsat);
			orb_epoch.mjd = md_julday(iyear,imon,iday);
			orb_epoch.sod = ih * 3600.0 + imin * 60.0 + dsec;
			if(lbdtime){ timinc(orb_epoch.mjd,orb_epoch.sod,14,&orb_epoch.mjd,&orb_epoch.sod);}
			for(int it = 0;it < nsat;it++){
				getline(in,line);
				sat_orb orb;
				sscanf(line.c_str(),"%s%d%lf%lf%lf%lf%lf%lf",orb.cprn,&orb.iode,orb.orb,orb.orb+1,orb.orb+2,orb.orb+3,
						orb.orb+4,orb.orb+5);
				orb_epoch.ORB[orb.cprn] = orb;
			}
			ORB_ALL.push_back(orb_epoch);
		}
		if(strstr(line.c_str(),"> CLOCK") != nullptr){
			clk clk_epoch;
			sscanf(line.c_str(),"%*s%*s%d%d%d%d%d%lf%d%d",&iyear,&imon,&iday,&ih,&imin,&dsec,&iv,&nsat);
			clk_epoch.mjd = md_julday(iyear,imon,iday);
			clk_epoch.sod = ih * 3600.0 + imin * 60.0 + dsec;
			if(lbdtime){ timinc(clk_epoch.mjd,clk_epoch.sod,14,&clk_epoch.mjd,&clk_epoch.sod);}
			for(int it = 0;it < nsat;it++){
				getline(in,line);
				sat_clk clk;
				sscanf(line.c_str(),"%s%d%lf%lf%lf",clk.cprn,&clk.iode,clk.clk,clk.clk + 1,clk.clk + 2);
				clk_epoch.CLK[clk.cprn] = clk;
			}
			CLK_ALL.push_back(clk_epoch);
		}
	}
	in.close();
}

void B2b::m_process(){
	m_loadorbclk();
	nav.v_openRnxEph(toString(mjd0) + ":" + toString(sod0));
	m_genclock();
	m_geneph();
}
void B2b::m_geneph(){
	FILE* fp = NULL;
	char filename[256] = {0};
	int mjd,iyear,imon,iday,ih,im,gweek,gd;
	double dx[6],length,A[3],Cr[3],R[3],dxx[6],sec;
	double sod,dintv = 900,dt,dtmin,dtlimit =120,xsat_out[6],xsat[6];
	/******** generate the ephemeris here *****************/
	cout << "begin to recover ephemeris ..." << endl;
	mjd = mjd0;
	sod = sod0;
	sod = ((int)(sod / dintv)) * dintv; /**force it to the epoch**/
	mjd2date(mjd,sod,&iyear,&imon,&iday,&ih,&im,&sec);
	gpsweek(iyear,imon,iday,&gweek,&gd);
	sprintf(filename,"cov%04d%d.sp3",gweek,gd);
	if(!(fp = fopen(filename,"w"))){cout << "can't open file " << filename << " to write!" << endl;exit(1);}
	m_wtephheader(fp);
	while(true){
		mjd2date(mjd,sod,&iyear,&imon,&iday,&ih,&im,&sec);
		fprintf(fp,"*   %4d%3d%3d%3d%3d%12.8lf\n",iyear,imon,iday,ih,im,sec);
		for(auto& cprn : cprn_list){
			/*************get the minus ssr*****************/
			dtmin = 9e9;
			orbit* ssr_selected = nullptr;
			for(auto& ssr : ORB_ALL){
				dt = (ssr.mjd - mjd) * 86400.0 + ssr.sod - sod;
				if(dtmin > fabs(dt)){
					dtmin = fabs(dt);
					ssr_selected = &ssr;
				}
			}
			if(dtmin > dtlimit || ssr_selected == nullptr) continue;
			/*************request for the brdm *****************/
			if(ssr_selected->ORB.find(cprn) == ssr_selected->ORB.end()) continue;
			if(!(nav.v_readOrbit(cprn.c_str(),mjd,sod,xsat_out,&ssr_selected->ORB[cprn].iode))){
//				cout << "ssr do not match ephemeris for " << cprn << " " << mjd << " " << sod << endl;
				continue;
			}
			/*************update the eph here *****************/
			memcpy(dx,ssr_selected->ORB[cprn].orb,sizeof(double) * 6);
			unit_vector(3,xsat_out,R,&length);
			cross(xsat_out,xsat_out + 3,Cr);
			unit_vector(3, Cr, Cr, &length);
			cross(Cr, R, A);
			unit_vector(3, A, A, &length);

//			unit_vector(3, xsat_out + 3, A, &length);
//			cross(xsat_out, xsat_out + 3, Cr);
//			unit_vector(3, Cr, Cr, &length);
//			cross(A, Cr, R);
//			unit_vector(3, R, R, &length);

			memset(dxx, 0, sizeof(dxx));
			dxx[0] = R[0] * dx[0] + A[0] * dx[1] + Cr[0] * dx[2];
			dxx[1] = R[1] * dx[0] + A[1] * dx[1] + Cr[1] * dx[2];
			dxx[2] = R[2] * dx[0] + A[2] * dx[1] + Cr[2] * dx[2];

			// velocity
			dxx[3] = R[0] * dx[3] + A[0] * dx[4] + Cr[0] * dx[5];
			dxx[4] = R[1] * dx[3] + A[1] * dx[4] + Cr[1] * dx[5];
			dxx[5] = R[2] * dx[3] + A[2] * dx[4] + Cr[2] * dx[5];

			for (int i = 0; i < 6; i++) {
				xsat[i] = xsat_out[i] - dxx[i];
			}
			double clock_sp3 = 999.0;
			/**** clock part *****/
			{
				double xclk_out,c[3],dclk;
				dtmin = 9e9;
				clk *ssr_selected = nullptr;
				for (auto &ssr : CLK_ALL) {
					dt = (ssr.mjd - mjd) * 86400.0 + ssr.sod - sod;
					if (dtmin > fabs(dt)) {
						dtmin = fabs(dt);
						ssr_selected = &ssr;
					}
				}
				if (dtmin < dtlimit && ssr_selected != nullptr){
					/*************request for the brdm *****************/
					if (ssr_selected->CLK.find(cprn) != ssr_selected->CLK.end()){
						if ((nav.v_readClk(cprn.c_str(), mjd, sod, &xclk_out, &ssr_selected->CLK[cprn].iode))) {
							/*************update the eph here *****************/
							memcpy(c, ssr_selected->CLK[cprn].clk, sizeof(double) * 3);
							dclk = c[0] + dtmin * (c[1] + dtmin * c[2]);
							clock_sp3 = xclk_out - dclk / VEL_LIGHT;
						}
					}
				}
			}
			/**************write it into orbit file *********/
			fprintf(fp,"P%s%14.6lf%14.6lf%14.6lf%14.6lf\n",cprn.c_str(), xsat[0] / 1000.0 , xsat[1] / 1000.0,xsat[2] / 1000.0,clock_sp3 != 999.0 ? clock_sp3 * 1E6 : 999.0);
		}
		timinc(mjd,sod,dintv,&mjd,&sod);
		if( (mjd - mjd1) * 86400.0 + sod - sod1 >= 0.0)
			break;
	}
	fclose(fp);
}
void B2b::m_geneph_brd(){
	FILE* fp = NULL;
	char filename[256] = {0};
	int mjd,iyear,imon,iday,ih,im,gweek,gd;
	double dx[6],length,A[3],Cr[3],R[3],dxx[6],sec;
	double sod,dintv = 900,dt,dtmin,dtlimit =120,xsat_out[6],xsat[6];
	/******** generate the ephemeris here *****************/
	cout << "begin to recover ephemeris ..." << endl;
	mjd = mjd0;
	sod = sod0;
	sod = ((int)(sod / dintv)) * dintv; /**force it to the epoch**/
	mjd2date(mjd,sod,&iyear,&imon,&iday,&ih,&im,&sec);
	gpsweek(iyear,imon,iday,&gweek,&gd);
	sprintf(filename,"cov%04d%d.sp3",gweek,gd);
	if(!(fp = fopen(filename,"w"))){cout << "can't open file " << filename << " to write!" << endl;exit(1);}
	m_wtephheader(fp);
	while(true){
		mjd2date(mjd,sod,&iyear,&imon,&iday,&ih,&im,&sec);
		fprintf(fp,"*   %4d%3d%3d%3d%3d%12.8lf\n",iyear,imon,iday,ih,im,sec);
		for(auto& cprn : cprn_list){
			/*************get the minus ssr*****************/
			if(!(nav.v_readOrbit(cprn.c_str(),mjd,sod,xsat_out))){
				cout << "ssr do not match ephemeris for " << cprn << " " << mjd << " " << sod << endl;
				continue;
			}
			for (int i = 0; i < 6; i++) {
				xsat[i] = xsat_out[i];
			}
			/**************write it into orbit file *********/
			fprintf(fp,"P%s%14.6lf%14.6lf%14.6lf%14.6lf\n",cprn.c_str(), xsat[0] / 1000.0 , xsat[1] / 1000.0,xsat[2] / 1000.0,999.0);
		}
		timinc(mjd,sod,dintv,&mjd,&sod);
		if( (mjd - mjd1) * 86400.0 + sod - sod1 >= 0.0)
			break;
	}
	fclose(fp);
}
void B2b::m_genclock_brd(){
	FILE* fp = NULL;
	char filename[128] = {0};
	int mjd,iyear,imon,iday,ih,im,gweek,gd;
	double sod,dintv = 300,dt,dtmin,dtlimit =120,xclk,xclk_out,dclk,c[3],sec;
	/******** generate the ephemeris here *****************/
	cout << "begin to recover clock ..." << endl;
	mjd = mjd0;
	sod = sod0;
	sod = ((int)(sod / dintv)) * dintv; /**force it to the epoch**/
	mjd2date(mjd,sod,&iyear,&imon,&iday,&ih,&im,&sec);
	gpsweek(iyear,imon,iday,&gweek,&gd);
	sprintf(filename,"cov%04d%d.clk",gweek,gd);
	if(!(fp = fopen(filename,"w"))){cout << "can't open file " << filename << " to write!" << endl;exit(1);}
	m_wtclkheader(fp);
	while(true){
		mjd2date(mjd,sod,&iyear,&imon,&iday,&ih,&im,&sec);
		for (auto &cprn : cprn_list){
			if (!(nav.v_readClk(cprn.c_str(),mjd,sod,&xclk_out))){
				cout << "ssr do not match ephemeris for " << cprn << " " << mjd << " " << sod << endl;
				continue;
			}
			/*************update the eph here *****************/
			xclk = xclk_out;
			fprintf(fp,"AS %3s%6d%3d%3d%3d%3d%10.6lf%3d%22.12e\r\n", cprn.c_str(), iyear, imon, iday, ih, im, sec, 1,xclk);
		}
		timinc(mjd,sod,dintv,&mjd,&sod);
		if( (mjd - mjd1) * 86400.0 + sod - sod1 >= 0.0)
			break;
	}
	fclose(fp);
}

void B2b::m_genclock(){
	FILE* fp = NULL;
	char filename[128] = {0};
	int mjd,iyear,imon,iday,ih,im,gweek,gd;
	double sod,dintv = 30,dt,dtmin,dtlimit =120,xclk,xclk_out,dclk,c[3],sec;
	/******** generate the ephemeris here *****************/
	cout << "begin to recover clock ..." << endl;
	mjd = mjd0;
	sod = sod0;
	sod = ((int)(sod / dintv)) * dintv; /**force it to the epoch**/
	mjd2date(mjd,sod,&iyear,&imon,&iday,&ih,&im,&sec);
	gpsweek(iyear,imon,iday,&gweek,&gd);
	sprintf(filename,"cov%04d%d.clk",gweek,gd);
	if(!(fp = fopen(filename,"w"))){cout << "can't open file " << filename << " to write!" << endl;exit(1);}
	m_wtclkheader(fp);
	while(true){
		mjd2date(mjd,sod,&iyear,&imon,&iday,&ih,&im,&sec);
		for (auto &cprn : cprn_list) {
			/*************get the minus ssr*****************/
			dtmin = 9e9;
			clk *ssr_selected = nullptr;
			for (auto &ssr : CLK_ALL) {
				dt = (ssr.mjd - mjd) * 86400.0 + ssr.sod - sod;
				if (dtmin > fabs(dt)) {
					dtmin = fabs(dt);
					ssr_selected = &ssr;
				}
			}
			if (dtmin > dtlimit || ssr_selected == nullptr)
				continue;
			/*************request for the brdm *****************/
			if (ssr_selected->CLK.find(cprn) == ssr_selected->CLK.end())
				continue;
			if (!(nav.v_readClk(cprn.c_str(),mjd,sod,&xclk_out,&ssr_selected->CLK[cprn].iode))){
				cout << "ssr do not match ephemeris for " << cprn << " " << mjd << " " << sod << endl;
				continue;
			}
			/*************update the eph here *****************/
			memcpy(c,ssr_selected->CLK[cprn].clk,sizeof(double) * 3);
			dclk = c[0] + dtmin * (c[1] + dtmin * c[2]);
			xclk = xclk_out - dclk / VEL_LIGHT;
			fprintf(fp,"AS %3s%6d%3d%3d%3d%3d%10.6lf%3d%22.12e\r\n", cprn.c_str(), iyear, imon, iday, ih, im, sec, 1,xclk);
		}
		timinc(mjd,sod,dintv,&mjd,&sod);
		if( (mjd - mjd1) * 86400.0 + sod - sod1 >= 0.0)
			break;
	}
	fclose(fp);
}
void B2b::m_wtclkheader(FILE* fp){
	char line[] ={
		"	 3.00           C                                       RINEX VERSION / TYPE\n"
		"     5.00                                                   INTERVAL\n"
		"PANDA-EXTCLK        WHU                 20200414            PGM / RUN BY / DATE\n"
		"G                   IGS14_2097.ATX                          SYS / PCVS APPLIED\n"
		"R                   IGS14_2097.ATX                          SYS / PCVS APPLIED\n"
		"C                   IGS14_2097.ATX                          SYS / PCVS APPLIED\n"
		"E                   IGS14_2097.ATX                          SYS / PCVS APPLIED\n"
		"    2    AR    AS                                           # / TYPES OF DATA\n"
		"WHU WUHAN UNIVERSITY                                        ANALYSIS CENTER\n"
		"  115                                                       # OF SOLN SATS\n"
		"G01 G02 G03 G05 G06 G07 G08 G09 G10 G11 G12 G13 G14 G15 G16 PRN LIST\n"
		"G17 G18 G19 G20 G21 G22 G24 G25 G26 G27 G28 G29 G30 G31 G32 PRN LIST\n"
		"R01 R02 R03 R04 R05 R07 R08 R09 R11 R12 R13 R14 R15 R16 R17 PRN LIST\n"
		"R18 R19 R20 R21 R22 R23 R24 C01 C02 C03 C06 C07 C08 C09 C10 PRN LIST\n"
		"C11 C12 C13 C14 C16 C19 C20 C21 C22 C23 C24 C25 C26 C27 C28 PRN LIST\n"
		"C29 C30 C32 C33 C34 C35 C36 C37 C38 C39 C40 C41 C42 C43 C44 PRN LIST\n"
		"C45 C46 C59 E01 E02 E03 E04 E05 E07 E08 E09 E11 E12 E13 E15 PRN LIST\n"
		"E19 E21 E24 E25 E26 E27 E30 E31 E33 E36                     PRN LIST\n"
		"                                                            END OF HEADER\n"};

	fprintf(fp,"%s",line);
	fflush(fp);
}
void B2b::m_wtephheader(FILE* fp){
	char line[] {
			"+  115   G01G02G03G05G06G07G08G09G10G11G12G13G14G15G16G17G18\n"
			"+        G19G20G21G22G24G25G26G27G28G29G30G31G32R01R02R03R04\n"
			"+        R05R07R08R09R11R12R13R14R15R16R17R18R19R20R21R22R23\n"
			"+        R24C01C02C03C06C07C08C09C10C11C12C13C14C16C19C20C21\n"
			"+        C22C23C24C25C26C27C28C29C30C32C33C34C35C36C37C38C39\n"
			"+        C40C41C42C43C44C45C46C59E01E02E03E04E05E07E08E09E11\n"
			"+        E12E13E15E19E21E24E25E26E27E30E31E33E36  0  0  0  0\n"
			"++        10  9  9  7 10  8 10 10  9 10  7 10 10  8 10  8 10\n"
			"++         8 10 10  8 10 10 10 10 10  9 10  8 10 10 10 10 10\n"
			"++        10 10 10 10 10 10 11 11 11 10 11 11 10 11 11 11 11\n"
			"++        11  9  9 11 10  9  9 10 10 10  9 10 10 10 10 10 10\n"
			"++        10 10 10 10 10  9 10 10  9 10 10 10  9 10 10  9  9\n"
			"++         9  9  9  9  9  9  9  9  9  9  9  9  9  9  9  9  9\n"
			"++         9  9  9  9  9  9  9  9  9  9  9  9  9  0  0  0  0\n"
			"%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n"
			"%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n"
			"%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n"
			"%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n"
			"%i    0    0    0    0      0      0      0      0         0\n"
			"%i    0    0    0    0      0      0      0      0         0\n" };

	int iyear,imon,iday,ih,imin,gweek,gd;
	double dsec,gsow;
	mjd2date(mjd0,sod0,&iyear,&imon,&iday,&ih,&imin,&dsec);
	gpsweek(iyear,imon,iday,&gweek,&gd);
	mjd2wksow(mjd0,sod0,&gweek,&gsow);
	fprintf(fp,"#dP%04d%3d%3d%3d%3d%12.8lf      96   u+U IGS14 FIT  WHU\n",iyear,imon,iday,ih,imin,dsec);
	fprintf(fp,"## %4d%16.8lf%15.8lf %d%16.8lf\n",gweek,gsow,30.0,mjd0,sod0);
	fprintf(fp,"%s",line);
	fflush(fp);
}

