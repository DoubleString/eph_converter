/*
 * Com.h
 *
 *  Created on: 2020年5月11日
 *      Author: dongdong
 */

#ifndef INCLUDE_COM_H_
#define INCLUDE_COM_H_
#include <iostream>
#include <vector>
#include "Const.h"
#include "Broadcast.h"

using namespace std;

extern int md_julday(int,int,int);

extern void xyzblh(double* x, double scale, double a0, double b0, double dx, double dy,
		double dz, double* geod);

void blhxyz(double* geod,double a0,double b0,double* x);

extern void rot_enu2xyz(double lat, double lon, double (*rotmat)[3]);

extern void m_wtclkhd(const char* filename,int mjd,double intv);

void brdtime(char* cprn,int *mjd,double *sod);

extern void matmpy(double* A, double* B, double* C, int row, int colA, int colB);

void wksow2mjd(int week,double sow,int* mjd,double* sod);

void gpsweek(int year, int month, int day, int* week, int* wd);

void mjd2doy(int jd, int* iyear, int* idoy);

extern int index_string(const char* src, char key);

extern int len_trim(const char* pStr);

void yr2year(int& yr);

void mjd2wksow(int mjd, double sod, int *week, double *sow);

double timdif(int jd2, double sod2, int jd1, double sod1);

extern void filleph(char* line,double ver);

char* substringEx(char* dest, const char* src, int start, int length);

int genAode(char csys,int mjd,double sod,double toe,int inade,GPSEPH* eph);

extern void timinc(int jd, double sec, double delt, int* jd1, double* sec1);

unsigned int getbdsiode(GPSEPH& bdseph);

unsigned long CRC24(long size, const unsigned char *buf) ;

void brdtime(char* cprn,int *mjd,double *sod);

void mjd2date(int jd, double sod, int* iyear, int* imonth, int* iday, int* ih,
		int* imin, double* sec);

void yeardoy2monthday(int iyear, int idoy, int* imonth, int* iday);

void unit_vector(int n, double* v, double* u, double* length);

void unit_vector(int n, double* v, double* u, double* length);

void cross(double* v1, double*v2, double* vout);


//extern void eph2pos(gtime_t time, const GPSEPH *eph, double *rs, double *dts, double *var);

#endif /* INCLUDE_COM_H_ */
