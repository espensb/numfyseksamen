#!/bin/bash
cd `dirname $0`
gcc=gcc
if command -v icc 
then
	gcc=icc
fi
$gcc
$gcc isingspinn.c -o isingspinn.o
time=`date +%s`
./isingspinn.o a > adata.dat
let time=`date +%s`-$time
echo "#Done in $time sec" >> adata.dat
echo "Done in $time sec" 
time=`date +%s`
./isingspinn.o b > bdata.dat
let time=`date +%s`-$time
echo "#Done in $time sec" >> bdata.dat
echo "Done in $time sec" 
