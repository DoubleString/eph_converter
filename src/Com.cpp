/*
 * Com.cpp
 *
 *  Created on: 2020年5月11日
 *      Author: dongdong
 */

#include <cmath>
#include "../include/Const.h"
#include "../include/Com.h"

using namespace std;

int md_julday(int iyear,int imonth,int iday){
	int iyr, result;
	int doy_of_month[12] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304,
			334 };
	if (iyear < 0 || imonth < 0 || iday < 0 || imonth > 12 || iday > 366
			|| (imonth != 0 && iday > 31)) {
		printf("iyear = %d,imonth = %d,iday = %d,incorrect argument",iyear,imonth,iday);
		printf("iyear = %d,imonth = %d,iday = %d,incorrect argument",iyear,imonth,iday);
		exit(1);
	}
	iyr = iyear;
	if (imonth <= 2)
		iyr -= 1;
	result = 365 * iyear - 678941 + iyr / 4 - iyr / 100 + iyr / 400 + iday;
	if (imonth != 0)
		result = result + doy_of_month[imonth - 1];
	return result;
}

void gpsweek(int year, int month, int day, int* week, int* wd) {
	int mjd;
	if (year != 0) {
		if (month != 0)
			mjd = md_julday(year, month, day);
		else
			mjd = md_julday(year, 1, 1) + day - 1;
	} else
		mjd = day;
	*week = (mjd - 44244) / 7;
	*wd = mjd - 44244 - (*week) * 7;
}

void mjd2doy(int jd, int* iyear, int* idoy) {
	*iyear = (jd + 678940) / 365;
	*idoy = jd - md_julday(*iyear, 1, 1);
	while (*idoy <= 0) {
		(*iyear)--;
		*idoy = jd - md_julday(*iyear, 1, 1) + 1;
	}
}

void wksow2mjd(int week,double sow,int* mjd,double* sod){
	if(mjd!=NULL)
		*mjd = (int)(sow/86400.0) + week*7 + 44244;
	if(sod!=NULL)
		*sod = fmod(sow,86400.0);
}

void blhxyz(double* geod,double a0,double b0,double* x){
	double a,b,e2,N,W;
	if(a0 == 0 || b0 == 0){
		a = 6378137.0;
		b = 298.257223563;  //alpha = (a-b)/a
	}else{
		a = a0;
		b = b0;
	}
	if(b <= 6000000)
		b = a - a / b;
	e2 = (a * a - b * b) / (a * a);

	W = sqrt(1 - e2 * pow(sin(geod[0]),2));
	N = a / W;
	x[0] = (N + geod[2]) * cos(geod[0]) * cos(geod[1]);
	x[1] = (N + geod[2]) * cos(geod[0]) * sin(geod[1]);
	x[2] = (N * (1-e2) + geod[2]) *sin(geod[0]);
}

void xyzblh(double* x, double scale, double a0, double b0, double dx, double dy,
		double dz, double* geod) {
	double a, b;
	long double xp, yp, zp, s, e2, rhd, rbd, n, zps, tmp1, tmp2;
	int i;

	if (a0 == 0 || b0 == 0) {
		a = 6378137.0;
		b = 298.257223563;  //alpha = (a-b)/a
	} else {
		a = a0;
		b = b0;
	}
	for (i = 0; i < 3; i++)
		geod[i] = 0;
	xp = x[0] * scale + dx;
	yp = x[1] * scale + dy;
	zp = x[2] * scale + dz;

	if (b <= 6000000)
		b = a - a / b;
	e2 = (a * a - b * b) / (a * a);
	s = sqrt(xp * xp + yp * yp);
	geod[1] = atan(yp / xp);
	if (geod[1] < 0) {
		if (yp > 0)
			geod[1] += PI;
		if (yp < 0)
			geod[1] += 2 * PI;
	} else {
		if (yp < 0)
			geod[1] += PI;
	}
	zps = zp / s;
	geod[2] = sqrt(xp * xp + yp * yp + zp * zp) - a;
	geod[0] = atan(zps / (1.0 - e2 * a / (a + geod[2])));
	n = 1;
	rhd = rbd = 1;
	while (rbd * n > 1e-4 || rhd > 1e-4) {
		n = a / sqrt(1.0 - e2 * sin(geod[0]) * sin(geod[0]));
		tmp1 = geod[0];
		tmp2 = geod[2];
		geod[2] = s / cos(geod[0]) - n;
		geod[0] = atan(zps / (1.0 - e2 * n / (n + geod[2])));
		rbd = ABS(tmp1 - geod[0]);
		rhd = ABS(tmp2 - geod[2]);
	}
}

void rot_enu2xyz(double lat, double lon, double (*rotmat)[3]) {
	double coslat, sinlat, coslon, sinlon;
	coslat = cos(lat - PI / 2);
	sinlat = sin(lat - PI / 2);
	coslon = cos(-PI / 2 - lon);
	sinlon = sin(-PI / 2 - lon);

	rotmat[0][0] = coslon;
	rotmat[0][1] = sinlon * coslat;
	rotmat[0][2] = sinlon * sinlat;
	rotmat[1][0] = -sinlon;
	rotmat[1][1] = coslon * coslat;
	rotmat[1][2] = coslon * sinlat;
	rotmat[2][0] = 0;
	rotmat[2][1] = -sinlat;
	rotmat[2][2] = coslat;
}

void matmpy(double* A, double* B, double* C, int row, int colA, int colB) {
	int i, j, k;
	double value;
	double* dest = (double*)calloc(row * colB,sizeof(double));
	for (i = 0; i < row; i++) {
		for (j = 0; j < colB; j++) {
			value = 0;
			for (k = 0; k < colA; k++){
				value += A[i * colA + k] * B[k* colB + j];
			}
			dest[i * colB + j] = value;
		}
	}
	for(i = 0;i <row;i++){
		for(j = 0;j < colB;j++){
			C[i * colB + j] = dest[i * colB + j];
		}
	}
	free(dest);
}

int index_string(const char* src, char key) {
	int len = strlen(src);
	int i;
	for (i = 0; i < len; i++) {
		if (src[i] == key)
			break;
	}
	if (i == len)
		return -1;
	else
		return i;
}

int len_trim(const char* pStr) {
	int length = strlen(pStr);
	int count = length;
	int i;
	for (i = length - 1; i >= 0; i--) {
		if (pStr[i] == '\0' || pStr[i] == '\n' || pStr[i] == '\r'
				|| pStr[i] == ' ')
			count--;
		else
			break;
	}
	return count;
}

void yr2year(int& yr) {
	if (yr > 1900)
		return;
	if (yr <= 30)
		yr += 2000;
	else
		yr += 1900;
}

void mjd2wksow(int mjd, double sod, int *week, double *sow) {
	*week = (int) ((mjd + sod / 86400.0 - 44244.0) / 7.0);
	*sow = (mjd - 44244.0 - *week * 7) * 86400.0 + sod;
}

double timdif(int jd2, double sod2, int jd1, double sod1) {
	return 86400.0 * (jd2 - jd1) + sod2 - sod1;
}

void filleph(char* line,double ver){
	int i,len,cpre = 4;
	char tmp[128];
	if(ver < 3.0)
		cpre = 3;
	len = len_trim(line);
	for(i = len;i < cpre + 19 * 4;i++){
		line[i] = ' ';
	}
	line[cpre + 19 * 4] = '\0';
	for (i = 0; i < 4; i++) {
		substringEx(tmp, line + cpre + 19 * i, 0, 19);
		if (len_trim(tmp) == 0) {
			line[cpre + 19 * i + 1] = '0';
		}
	}
}

// start: started with index of zero
char* substringEx(char* dest, const char* src, int start, int length) {
	int i, j = 0;
	int len = strlen(src);
	if (start < 0 || start >= len || start + length > len) {
		dest[0] = '\0';
		return NULL;
	}
	for (i = start; i < start + length; i++) {
		dest[j] = src[i];
		j++;
	}
	dest[j] = '\0';
	return dest;
}

int genAode(char csys,int mjd,double sod,double toe,int inade,GPSEPH* eph){
	int mjd_r,ret = -1;
	double sod_r;
	if(csys == 'G' || csys == 'J')
		ret = inade;
	else if(csys == 'R'){
		timinc(mjd,sod, 19.0 - 37,&mjd_r,&sod_r);
		//ret = NINT(sod / 900.0) + 1;
		//IGMAS:
		ret = NINT(fmod(sod_r + 10800.0, 86400.0) / 900.0);
	}
	else if(csys == 'E'){
		// ret = NINT(sod / 600.0) + 1;
		//IGMAS:
		ret = inade;
	}
	else if(csys == 'C'){
		// ret = NINT(sod / 1800.0) + 1;
		//IGMAS:
		//ret = NINT(fmod(toe, 86400.0) / 450.0);
//		ret = getbdsiode(*eph);
		ret = (int)eph->iodc;
	}
	return ret;
}

void timinc(int jd, double sec, double delt, int* jd1, double* sec1) {
	*sec1 = sec + delt;
	int inc = (int) (*sec1 / 86400.0);
	*jd1 = jd + inc;
	*sec1 = *sec1 - inc * 86400.0;
	if (*sec1 >= 0)
		return;
	*jd1 = *jd1 - 1;
	*sec1 = *sec1 + 86400;
}

unsigned int getbdsiode(GPSEPH& bdseph) {
	unsigned char buffer[80];
	int size = 0;
	int numbits = 0;
	long long bitbuffer = 0;
	unsigned char *startbuffer = buffer;

	BDSADDBITSFLOAT(14, bdseph.idot, M_PI/(double)(1<<30)/(double)(1<<13));
	BDSADDBITSFLOAT(11, bdseph.a2,
			1.0 / (double )(1 << 30) / (double )(1 << 30) / (double )(1 << 6));
	BDSADDBITSFLOAT(22, bdseph.a1,
			1.0 / (double )(1 << 30) / (double )(1 << 20));
	BDSADDBITSFLOAT(24, bdseph.a0, 1.0 / (double )(1 << 30) / (double )(1 << 3));
	BDSADDBITSFLOAT(18, bdseph.crs, 1.0 / (double )(1 << 6));
	BDSADDBITSFLOAT(16, bdseph.dn, M_PI/(double)(1<<30)/(double)(1<<13));
	BDSADDBITSFLOAT(32, bdseph.m0, M_PI/(double)(1<<30)/(double)(1<<1));
	BDSADDBITSFLOAT(18, bdseph.cuc,
			1.0 / (double )(1 << 30) / (double )(1 << 1));
	BDSADDBITSFLOAT(32, bdseph.e, 1.0 / (double )(1 << 30) / (double )(1 << 3));
	BDSADDBITSFLOAT(18, bdseph.cus,
			1.0 / (double )(1 << 30) / (double )(1 << 1));
	BDSADDBITSFLOAT(32, bdseph.roota, 1.0 / (double )(1 << 19));
	BDSADDBITSFLOAT(18, bdseph.cic,
			1.0 / (double )(1 << 30) / (double )(1 << 1));
	BDSADDBITSFLOAT(32, bdseph.omega0, M_PI/(double)(1<<30)/(double)(1<<1));
	BDSADDBITSFLOAT(18, bdseph.cis,
			1.0 / (double )(1 << 30) / (double )(1 << 1));
	BDSADDBITSFLOAT(32, bdseph.i0, M_PI/(double)(1<<30)/(double)(1<<1));
	BDSADDBITSFLOAT(18, bdseph.crc, 1.0 / (double )(1 << 6));
	BDSADDBITSFLOAT(32, bdseph.omega, M_PI/(double)(1<<30)/(double)(1<<1));
	BDSADDBITSFLOAT(24, bdseph.omegadot, M_PI/(double)(1<<30)/(double)(1<<13));
	BDSADDBITS(5, 0); // the last byte is filled by 0-bits to obtain a length of an integer multiple of 8

	return CRC24(size, startbuffer);
}

unsigned long CRC24(long size, const unsigned char *buf) {
	unsigned long crc = 0;
	int ii;
	while (size--) {
		crc ^= (*buf++) << (16);
		for (ii = 0; ii < 8; ii++) {
			crc <<= 1;
			if (crc & 0x1000000) {
				crc ^= 0x01864cfb;
			}
		}
	}
	return crc;
}

void brdtime(char* cprn,int *mjd,double *sod){
	switch (cprn[0]) {
	case 'C':
	case 'B':
		timinc(*mjd, *sod, 14.0, mjd, sod);
		break;
	case 'R':
		//timinc(*mjd, *sod, Taiutc::s_getInstance()->m_getTaiutc(*mjd) - 19.0, mjd, sod);
		break;
	}
}

/* broadcast ephemeris to satellite position and clock bias --------------------
* compute satellite position and clock bias with broadcast ephemeris (gps,
* galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          eph_t *eph       I   broadcast ephemeris
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [1],[7],[8]
*          satellite clock includes relativity correction without code bias
*          (tgd or bgd)
*-----------------------------------------------------------------------------*/
/*extern void eph2pos(gtime_t time, const GPSEPH *eph, double *rs, double *dts,
                    double *var)
{
	static double A_ref_igso=42162200, A_ref_MEO=27906100;
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
	double na,delta_na,A0,Ak,n0;
    int n,sys,prn;

    //trace(4,"eph2pos : time=%s sat=%2d\n",time_str(time,3),eph->sat);

    if (eph->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
        return;
    }
    tk=timediff(time,eph->toe);

    switch ((sys=satsys(eph->sat,&prn))) {
        case SYS_GAL: mu=MU_GAL; omge=OMGE_GAL; break;
        case SYS_CMP: mu=MU_CMP; omge=OMGE_CMP; break;
        default:      mu=MU_GPS; omge=OMGE;     break;
    }
	if(sys==SYS_CMP&&eph->signal_idx>6){
		A0=((prn>=38&&prn<=40||(prn>=59&&prn<=63))?A_ref_igso:A_ref_MEO)+eph->delta_A;
		Ak=A0+eph->A_DOT*tk;
		//eph->A=Ak;
		n0=sqrt(mu/A0/A0/A0);
		delta_na=eph->deln+0.5*eph->delta_n_dot*tk;
		na=n0+delta_na;
		M=eph->M0+na*tk;
	}
	else{
		M=eph->M0+(sqrt(mu/(eph->A*eph->A*eph->A))+eph->deln)*tk;
	}

    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<100;n++) {
        Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
    }
	if(n>=100){rs[0]=rs[1]=rs[2]=0.0;*dts=0.0;return;}
    sinE=sin(E); cosE=cos(E);

    //trace(4,"kepler: sat=%2d e=%8.5f n=%2d del=%10.3e\n",eph->sat,eph->e,n,E-Ek);

    u=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;
    r=((sys==SYS_CMP&&eph->signal_idx>6)?Ak:eph->A)*(1.0-eph->e*cosE);
    i=eph->i0+eph->idot*tk;
    sin2u=sin(2.0*u); cos2u=cos(2.0*u);
    u+=eph->cus*sin2u+eph->cuc*cos2u;
    r+=eph->crs*sin2u+eph->crc*cos2u;
    i+=eph->cis*sin2u+eph->cic*cos2u;
    x=r*cos(u); y=r*sin(u); cosi=cos(i);

    if (sys==SYS_CMP&&(prn<=5||prn>=59&&prn<=63)) {
        O=eph->OMG0+eph->OMGd*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        xg=x*cosO-y*cosi*sinO;
        yg=x*sinO+y*cosi*cosO;
        zg=y*sin(i);
        sino=sin(omge*tk); coso=cos(omge*tk);
        rs[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
        rs[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
        rs[2]=-yg*SIN_5+zg*COS_5;
    }
    else {
        O=eph->OMG0+(eph->OMGd-omge)*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        rs[0]=x*cosO-y*cosi*sinO;
        rs[1]=x*sinO+y*cosi*cosO;
        rs[2]=y*sin(i);
    }
    tk=timediff(time,eph->toc);
    *dts=eph->f0+eph->f1*tk+eph->f2*tk*tk;

	*dts -= 2.0*sqrt(mu*((sys == SYS_CMP&&eph->signal_idx>6) ? Ak : eph->A))*eph->e*sinE / SQR(CLIGHT);

    *var=var_uraeph(eph->sva);
}*/

void mjd2date(int jd, double sod, int* iyear, int* imonth, int* iday, int* ih,
		int* imin, double* sec) {
	int doy = 0;
	mjd2doy(jd, iyear, &doy);
	yeardoy2monthday(*iyear, doy, imonth, iday);

	*ih =  static_cast<int>( sod / 3600.0);
	*imin = static_cast<int> ((sod - (*ih) * 3600.0) / 60.0);
	*sec = sod - (*ih) * 3600.0 - (*imin) * 60.0;
}

void yeardoy2monthday(int iyear, int idoy, int* imonth, int* iday) {
	int days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	int id, i;
	if ((iyear % 4 == 0 && iyear % 100 != 0) || iyear % 400 == 0)
		days_in_month[1] = 29;
	id = idoy;
	for (i = 0; i < 12; i++) {
		id = id - days_in_month[i];
		if (id > 0)
			continue;
		*iday = id + days_in_month[i];
		*imonth = i + 1;
		break;
	}
}

void unit_vector(int n, double* v, double* u, double* length) {

	int i;
	*length = 0.0;
	for (i = 0; i < n; i++) {
		*length = *length + v[i] * v[i];
	}
	(*length) = sqrt(*length);
	for (i = 0; i < n; i++)
		u[i] = v[i] / (*length);
}

void cross(double* v1, double*v2, double* vout) {
	double tmp[3];
	tmp[0] = v1[1] * v2[2] - v1[2] * v2[1];
	tmp[1] = v1[2] * v2[0] - v1[0] * v2[2];
	tmp[2] = v1[0] * v2[1] - v1[1] * v2[0];

	vout[0] = tmp[0];
	vout[1] = tmp[1];
	vout[2] = tmp[2];
}

