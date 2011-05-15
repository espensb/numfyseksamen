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
./isingspinn.o a > data.dat
let time=`date +%s`-$time
echo "Done in $time sec\n" 
gnuplot doplot.gnu
