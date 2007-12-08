#!/usr/bin/perl


foreach $f ( @ARGV ) {
    open(FH,"$f");
    while(<FH>) {
	s/\s+//g;
	if ( ! /^\s*$/ ) {
	    $l .= uc($_).",";
	}
    }
    close(FH);
}
chop($l);
print $l;
