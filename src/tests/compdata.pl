#!/usr/bin/perl

# compare limit in percent
$limit = 0.1;

open(RH,"out.ref");
open(FH,"out.0.gpl");

$failed = 0;
while(<RH>) {
    $r = $_;
    $c = <FH>;
    if ( ! /^\s*\#/ && ! /^\s*$/ ) {
	$c =~ /(\S+)$/;
	$val = $1;
	$r =~ /(\S+)$/;
	$rval = $1;
	if ( $rval != 0. ) {
	    $diff = abs(($val - $rval)/$rval)* 100.;
	    if ( $diff > $limit ) { $failed = 1; last; }
	} else {
	   if ( $val != 0. ) {  $failed = 1; last; }	
	}

    } 
}


close(RH);
close(FH);

exit($failed);
