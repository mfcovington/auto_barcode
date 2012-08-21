#!/usr/bin/env perl
# basic_barcode_split.v01.pl
# Mike Covington
# created: 2012-02-21
#
# Description: 
#
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Printer;

#options/defaults
my ( $fq_in, $barcode, $sample_id, $list, $help );
my $out_dir = "./";
my $options = GetOptions(
    "fq_in=s"     => \$fq_in,
    "barcode=s"   => \$barcode,
    "sample_id=s" => \$sample_id,
    "list"        => \$list,
    "out_dir=s"   => \$out_dir,
    "help"        => \$help,
);

#help/usage
my $prog = basename($0);
print_usage() and exit if $help;
print_usage() and exit unless defined $fq_in and defined $barcode;
print_usage() and exit unless defined $sample_id or defined $list;

#gather barcodes to search for
my @barcode_list;
if ($list) {
    open my $barcode_list_fh, '<', $barcode;
    while ( my $line = <$barcode_list_fh> ) {
        chomp $line;
        push @barcode_list, [ split /\t/, $line ];
    }
    close $barcode_list_fh;
    # $barcode = basename($barcode);
}
else { push @barcode_list, [ $sample_id, $barcode ]; }

my ($filename, $directories, $suffix) = fileparse($fq_in, "(.fq)|(.fastq)");
my $fq_out = $directories . $filename . ".$barcode" . $suffix;

open (FQ_IN, $fq_in) or die "Can't open $fq_in";
open (FQ_OUT, ">$fq_out") or die "Can't open $fq_out";

my $count = 0;
while (my $read_id = <FQ_IN>) {
	my $seq = <FQ_IN>;
	if ($seq =~ m/^$barcode/) {
		print FQ_OUT $read_id . $seq . <FQ_IN> . <FQ_IN>;
		$count++;
	}else{
		<FQ_IN>, <FQ_IN>;
	}
}

print "\n\tTotal # of reads with $barcode barcode = ", $count, "\n\n";

close (FQ_IN);
close (FQ_OUT);



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
  -o, --out_dir   DIR       Output file is saved in the specified directory
                              (or same directory as IN.FASTQ, if --out_dir is not used)

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


###CHALLENGES:
# 1. Alter the script to trim the barcode from the reads after parsing. (hint: don't forget to trim the quality scores and the 'T' following the barcode)
#		Hint: The substr function is handy sometimes.
# 2. Alter the script so that it reads multiple barcodes from a file (barcode.csv) or the command line and parses for each barcode.
#		Hint #1: Use the split function to parse the lines of the barcode file.
#		Hint #2: Of the two main ways to organize your program, one should be much faster than the other (especially for large FASTQ files).
#				 Think about what might be the best approach before you start writing your script. 