#!/usr/bin/perl
open($adata,"adata.dat");
$write = 0;
while (<$adata>) {
	if (/#N = (\d+)/) {
		close $out if ($out);
		$N = $1;
		open $out, ">data/adata$N.dat";
		next;
	}
	print $out $_;
}
close $out;
