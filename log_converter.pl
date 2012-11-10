#!/usr/bin/env perl
# log_converter.pl
# Mike Covington
# created: 2012-10-07
#
# Description: 
#
use strict;
use warnings;
use autodie;
use feature 'say';

die unless scalar @ARGV == 1;
my $log_file = $ARGV[0];
my $out_file = "$log_file.tsv";
open my $log_fh, "<", $log_file;
open my $out_fh, ">", $out_file;

say $out_fh "barcode\tcount\tmatched";
<$log_fh>;
while (<$log_fh>) {
    my ( $barcode, $count, $match ) =
      $_   =~ m/ ([ACTGN]+) \s+ ([\d,]+) .+ % \s* ([^\s]*) /ix;
    $count =~ s/,//g;
    if ( length $match > 0 ) {
        $match = "matched";
    }
    else { $match = "unmatched" }
    say $out_fh join "\t", $barcode, $count, $match;
}
close $log_fh;
close $out_fh;