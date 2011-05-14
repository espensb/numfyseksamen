#!/usr/bin/perl
use File::Basename;
chdir dirname($0);

`gcc isingspinn.c -o isingspinn.o`;
`./isingspinn.o > data.dat`;
`gnuplot doplot.gnu`;
`open plot.png`;
