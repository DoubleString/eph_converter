/*
 * main.cpp
 *
 *  Created on: 2020年5月11日
 *      Author: dongdong
 */

#include "../include/Const.h"
#include "../include/Com.h"
#include "../include/Broadcast.h"
#include "../include/B2bDecoder.h"

using namespace std;

void printUsages(){
	printf("usage:BDS-B2b-Decoder year month day hour minute second seslen ssr-file navfile -bdtime -help\n");
}
int main(int argc,char* args[]){
	B2b b2b;
	if(argc < 10){
		printUsages();
		exit(1);
	}
	b2b.lbdtime = false;
	for(int iargc= 0;iargc < argc;iargc++){
		if(strstr(args[iargc],"-bdtime"))
			b2b.lbdtime = true;
		if(strstr(args[iargc],"-help")){
			printUsages();
			exit(1);
		}
	}
	int iyear,imon,iday,ihour,imin;
	double dsec,dseslen;
	string obsfile,ephfile;
	iyear=atoi(args[1]);
	imon=atoi(args[2]);
	iday=atoi(args[3]);
	ihour=atoi(args[4]);
	imin=atoi(args[5]);
	dsec=atof(args[6]);
	dseslen = atof(args[7]);
	obsfile=string(args[8]);
	b2b.mjd0 = md_julday(iyear,imon,iday);
	b2b.sod0 = ihour * 3600.0 + imin * 60.0 + dsec;
	timinc(b2b.mjd0,b2b.sod0,dseslen,&b2b.mjd1,&b2b.sod1);
	b2b.nav.navfilename=string(args[9]);
	strcpy(b2b.File,obsfile.c_str());

	b2b.m_process();
	return 0;
}
