#!/usr/bin/perl
use Getopt::Std;

%options=();
getopts("s",\%options);
if ( defined ( $options{s} ) ) {
    $sep = " ";
} else {
    $sep = ",";
}
foreach $f ( @ARGV ) {
    open(FH,"$f");
    while(<FH>) {
	s/\s+//g;
	if ( ! /^\s*$/ ) {
	    if ( -f $_.".f90" ) {
		if ($sep eq ",") {
		    $l .= uc($_).$sep;
		} else {
		    $l .= lc($_).$sep;
		}
	    }
	}
    }
    close(FH);
}
chop($l);
print $l;
