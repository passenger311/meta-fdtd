#!/usr/bin/perl

    sub getnum {
        use POSIX qw(strtod);
        my $str = shift;
        $str =~ s/^\s+//;
        $str =~ s/\s+$//;
        $! = 0;
        my($num, $unparsed) = strtod($str);
        if (($str eq '') || ($unparsed != 0) || $!) {
            return undef;
        } else {
            return $num;
        } 
    } 

    sub is_numeric { defined &getnum } 

# compare limit in percent
$limit = 0.1;

open(RH,"out.ref");
open(FH,"out.0.gpl");

$failed = 0;
while(<RH>) {
    $r = $_;
    if ( $r =~ /^\s*$/ || $r =~ /^\s*\#/ ) { next; } 
    $c = <FH>;
    while ( $c =~ /^\s*$/ || $c =~ /^\s*\#/ ) {
	$c = <FH>;
    }
    if ( ! /^\s*$/ && ! /^\#/ ) {
	$c =~ /(\S+)$/;
	$val = $1;
	$r =~ /(\S+)$/;
	$rval = $1;
	if (! is_numeric $val ) { $failed = 1; last; } 
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
