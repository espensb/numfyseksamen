#!/usr/bin/perl
use File::Basename;
chdir dirname($0);

`gcc isingspinn.c -o isingspinn.o`;
$time = time;
`./isingspinn.o > data.dat`;
$time = time - $time;
print "Done in $time sec\n"; 
`gnuplot doplot.gnu`;
`open plot.png`;
