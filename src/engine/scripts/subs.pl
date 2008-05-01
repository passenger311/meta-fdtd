#!/usr/bin/perl
# -i bak



for ( $i =2; $i<=$#ARGV; $i++ ) {
$out = "$ARGV[$i].new";
open(GH,">$out");
open(FH,$ARGV[$i]);
while(<FH>) { 
    s/$ARGV[0]/$ARGV[1]/g; 
    print GH $_;

} 
close(FH);
close(GH);
`mv $ARGV[$i] $ARGV[$i].old`;
`mv $ARGV[$i].new $ARGV[$i]`;
}
