#!/bin/bash

main(){
	if [ $# -ne 5 ];then
		printUsage
		exit 1
	fi
	year=$1
	month=$2
	day=$3
	
	b2bfile=$4
	ephfile=$5
	
	./b2b-decoder -ssr $b2bfile -out ssr_decoder_b2b
	./ephconverter $year $month $day 0 0 0 86400 ssr_decoder_b2b $ephfile -bdtime
 
}
printUsage(){
	echo "-----------------------------------------------------"
	echo " GNSS B2B2SP3 SCRIPT FOR BAMBOO PROGRAM[CREATED BY JUNTAO AT WUHAN-UNIVERSITY]"
	echo " Usage: b2b2sp3.sh year month day b2bfile ephfile[rinex/18-pam for bd3]"
	echo "-----------------------------------------------------"
}
######################################################################
main "$@"