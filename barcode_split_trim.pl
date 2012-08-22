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
my ( $fq_in, $barcode, $id, $list, $outdir, $notrim, $autoprefix, $autosuffix, $help );
my $prefix = "";
my $suffix = "";
my $options = GetOptions(
    "fq_in=s"    => \$fq_in,
    "barcode=s"  => \$barcode,
    "id=s"       => \$id,
    "list"       => \$list,
    "outdir=s"   => \$outdir,
    "autoprefix" => \$autoprefix,
    "prefix=s"   => \$prefix,
    "suffix=s"   => \$suffix,
    "autosuffix" => \$autosuffix,
    "notrim"     => \$notrim,
    "help"       => \$help,
);

#help/usage
my $prog = basename($0);
print_usage() and exit if $help;
print_usage() and exit unless defined $fq_in and defined $barcode;
print_usage() and exit unless defined $id or defined $list;

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
else { $barcode_table{$barcode} = [ $id, 0 ]; }

#check lengths of barcodes
my $min_length = min map { length } keys %barcode_table;
my $max_length = max map { length } keys %barcode_table;
die "Unexpected variation in barcode length (min=$min_length, max=$max_length)"
  unless $min_length == $max_length;
my $barcode_length = $max_length;

#open all filehandles
my ( $filename, $directory, $filesuffix ) = fileparse( $fq_in, ".f(ast)?q" );
$prefix .= "." unless $prefix eq "";
$prefix = join ".", $filename, $prefix if $autoprefix;
$suffix = "." . $suffix unless $suffix eq "";
$directory = $outdir if defined $outdir;
open my $fq_in_fh, "<", $fq_in;
for ( keys %barcode_table ) {
    my $temp_suffix = $suffix;
    $temp_suffix = join ".", $suffix, $_ if $autosuffix;
    my $fq_out = $directory . $prefix . $barcode_table{$_}->{id} . $temp_suffix . ".fq";
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
open my $log_fh, ">", $directory . join ".", "barcode_splitting", "fq_" . $filename, "bar_" . $barcode, "log";
say $log_fh "Barcode spiltting sumary for $fq_in";
say $log_fh "-----------------------------" . "-" x length $fq_in;
say $log_fh "barcode\tid\tcount";
map {
    say $log_fh $_ . "\t"
      . $barcode_table{$_}->{id} . "\t"
      . commify( $barcode_table{$_}->{count} )
} keys %barcode_table;
say $log_fh "matched\t" . commify($total_matched);
say $log_fh "none\t"    . commify($total_unmatched);

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

OPTIONS
  -h, --help                 Print this help message
  -f, --fq_in     IN.FASTQ   Extract reads from specified file    
  -b, --barcode   BARCODE    Specify barcode or list of barcodes to extract
  -i, --id        SAMPLE_ID  Sample ID (not needed if using list of barcodes)
  -l, --list                 Indicates --barcode is a list of barcodes in a file
  -n, --notrim               Split without trimming barcodes off
  -o, --outdir    DIR        Output file is saved in the specified directory
                              (or same directory as IN.FASTQ, if --outdir is not used)

NAMING OPTIONS
  --autoprefix               Append FASTQ file name onto output
  --autosuffix               Append barcode onto output 
  -p, --prefix    PREFIX     Add custom prefix to output
  -s, --suffix    SUFFIX     Add custom suffix to output

OUTPUT
  An output file in fastq format is written for each barcode to the directory
  containing IN.FASTQ, unless an output directory is specified.
  The default name of the output file is SAMPLE_ID.fq. The output names can be
  customized using the Naming Options.

EXAMPLES
  $prog -f kitten_DNA.fq -b GACTG -i Charlotte
  $prog --fq_in kitten_DNA.fastq --barcode barcode.file --list
  $prog --help

EOF
}

exit;
