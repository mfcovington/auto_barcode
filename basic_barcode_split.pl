#!/usr/bin/env perl
# basic_barcode_split.pl
# Mike Covington
# created: 2012-02-21
#
# Description: 
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use File::Basename;
use List::Util qw(min max);

###TODO:
#check that barcodes are comprosed of ACGT
#update usage statement
#incorporate 'barcode_psychic.pl' functionality

#options/defaults
my ( $fq_in, $barcode, $sample_id, $list, $outdir, $notrim, $help );
my $options = GetOptions(
    "fq_in=s"     => \$fq_in,
    "barcode=s"   => \$barcode,
    "sample_id=s" => \$sample_id,
    "list"        => \$list,
    "outdir=s"    => \$outdir,
    "notrim"      => \$notrim,
    "help"        => \$help,
);

#help/usage
my $prog = basename($0);
print_usage() and exit if $help;
print_usage() and exit unless defined $fq_in and defined $barcode;
print_usage() and exit unless defined $sample_id or defined $list;

#gather barcodes to search for
my %barcode_table;
if ($list) {
    open my $barcode_list_fh, '<', $barcode;
    %barcode_table =
      map {
        chomp;
        my @delim = split /\t/;
        ( $delim[0], { 'id' => $delim[1], 'count' => 0 } )
      } <$barcode_list_fh>;
    close $barcode_list_fh;
}
else { $barcode_table{$barcode} = [ $sample_id, 0 ]; }

#check lengths of barcodes
my $min_length = min map { length } keys %barcode_table;
my $max_length = max map { length } keys %barcode_table;
die "Unexpected variation in barcode length (min=$min_length, max=$max_length)"
  unless $min_length == $max_length;
my $barcode_length = $max_length;

#open all filehandles
my ( $filename, $directory, $suffix ) = fileparse( $fq_in, ".f(ast)?q" );
$directory = $outdir if defined $outdir;
open my $fq_in_fh, "<", $fq_in;
for ( keys %barcode_table ) {
    my $fq_out = $directory . $barcode_table{$_}->{id} . ".fq";
    open $barcode_table{$_}->{fh}, ">", $fq_out;
}

#split and trim
my $total_matched   = 0;
my $total_unmatched = 0;
while ( my $read_id = <$fq_in_fh> ) {
    my $seq = <$fq_in_fh>;
    my $cur_barcode = substr $seq, 0, $barcode_length;
    if ( /^$cur_barcode/ ~~ %barcode_table ) {
        $seq = substr $seq, $barcode_length + 1 unless $notrim;
        my $qual_id = <$fq_in_fh>;
        my $qual    = <$fq_in_fh>;
        $qual = substr $qual, $barcode_length + 1 unless $notrim;
        print { $barcode_table{$cur_barcode}->{fh} }
          $read_id . $seq . $qual_id . $qual;
        $barcode_table{$cur_barcode}->{count}++;
        $total_matched++;
    }
    else {
        <$fq_in_fh>, <$fq_in_fh>;
        $total_unmatched++;
    }
}
map { close $barcode_table{$_}->{fh} } keys %barcode_table;

#summary
say "Barcode spiltting sumary for $fq_in";
say "-----------------------------" . "-" x length $fq_in;
say "barcode\tsample_id\tcount";
map {
    say $_ . "\t"
      . $barcode_table{$_}->{id} . "\t"
      . commify( $barcode_table{$_}->{count} )
} keys %barcode_table;
say "matched\t" . commify($total_matched);
say "none\t"    . commify($total_unmatched);

exit;

sub commify {
    local $_  = shift;
    1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
    return $_;
}

sub print_usage {
    warn <<"EOF";

USAGE
  $prog [options] -f IN.FASTQ -b BARCODE [-o OUTPUT_DIR]

DESCRIPTION
  Extracts fastq reads for specified barcode(s)

OPTIONS  ###THIS NEEDS TO BE UPDATED###
  -h, --help                Print this help message
  -f, --fq_in     IN.FASTQ  Extract reads from specified file    
  -b, --barcode   BARCODE   Specify barcode or list of barcodes to extract
  -s, --sample_id
  -l, --list                Indicates --barcode is a list of barcodes in a file
  -o, --outdir   DIR       Output file is saved in the specified directory
                              (or same directory as IN.FASTQ, if --outdir is not used)

OUTPUT  ###THIS NEEDS TO BE UPDATED###
  An output file in fastq format is written to the current directory, 
  unless an output directory is specified.
  The name of the output file is BARCODE.fastq.

EXAMPLES  ###THIS NEEDS TO BE UPDATED###
  $prog -f ITAG2.3_cds_SHORTNAMES.fastq -g Solyc10g044670 -o seq_directory
  $prog --fq_in ITAG2.3_cds_SHORTNAMES.fastq --barcode barcode.file --list
  $prog --help

EOF
}

exit;
