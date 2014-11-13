#!/usr/bin/env perl
# barcode_split_trim.pl
# Mike Covington (Maloof Lab, UC-Davis)
# https://github.com/mfcovington/auto_barcode
#
# Description:
# - Extracts fastq reads for specified barcode(s) from one or multiple FASTQ files
# - Writes helpful logs with barcode stats
# - Plots an exquisite summary of matched and unmatched barcodes
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use File::Basename;
use File::Path 'make_path';
use List::MoreUtils 'part';
use List::Util qw(min max);
use Statistics::Descriptive;
use Statistics::R;
use Text::Levenshtein::XS 'distance';
use Text::Table;

# TODO: incorporate more 'barcode_psychic.pl' functionality (warnings/suggestions)
# TODO: fuzzy matching
# TODO: Change 'matched/unmatched' to 'expected/unexpected'

my $current_version = "v2.0.0";

#options/defaults
my $mismatches_ok = 0;
my ($barcode,     $id,     $list,  $outdir,
    $indexed, $notrim, $stats, $autoprefix,
    $autosuffix,  $help,   $version
);
my $prefix  = "";
my $suffix  = "";
my $options = GetOptions(
    "barcode=s"    => \$barcode,
    "id=s"         => \$id,
    "list"         => \$list,
    "mismatches=s" => \$mismatches_ok,
    "outdir=s"     => \$outdir,
    "indexed"      => \$indexed,
    "autoprefix"   => \$autoprefix,
    "prefix=s"     => \$prefix,
    "suffix=s"     => \$suffix,
    "autosuffix"   => \$autosuffix,
    "notrim"       => \$notrim,
    "stats"        => \$stats,
    "help"         => \$help,
    "version"      => \$version,
);
my @all_fq_files = grep { /f(ast)?q$/i } @ARGV;
my ( $fq_files, $fq_indexes )
    = $indexed ? part { my $i++ % 2 } @all_fq_files : \@all_fq_files;

validate_options( $version, $help, $barcode, $id, $list, $mismatches_ok,
    $fq_files, $fq_indexes, $indexed );

my $barcode_table = get_barcodes( $list, $barcode, $id );

my $barcode_length = validate_barcodes($barcode_table);

my ( $directory, $filename, $barcode_name )
    = parse_filenames( $fq_files, $outdir, $barcode );

make_path($directory);

my $unmatched_fh
    = open_fq_fhs( $barcode_table, $directory, $filename, $barcode_name,
                   $prefix, $suffix, $autoprefix, $autosuffix ) unless $stats;

my ( $total_matched, $total_unmatched, $barcodes_obs )
    = split_trim_barcodes( $fq_files, $fq_indexes, $indexed,
    $barcode_table, $barcode_length, $unmatched_fh, $notrim, $stats );

close_fq_fhs( $barcode_table, $unmatched_fh ) unless $stats;

my $total_count = $total_matched + $total_unmatched;

summarize_observed_barcodes(
    $barcode_table, $barcodes_obs, $total_count,
    $directory,     $filename,     $barcode_name
);

summarize_counts(
    $barcode_table, $fq_files, $total_count,
    $total_matched, $total_unmatched, $directory, $filename, $barcode_name
);

plot_summary( $barcodes_obs, $barcode_table, $directory, $id );

exit;

sub percent {
    local $_ = shift;
    return sprintf( "%.1f", $_ * 100 ) . "%";
}

sub round {
    local $_ = shift;
    return int( $_ + 0.5 );
}

sub commify {
    local $_ = shift;
    1 while s/^([-+]?\d+)(\d{3})/$1,$2/;
    return $_;
}

sub validate_options {
    my ($version,  $help,       $barcode,
        $id,       $list,       $mismatches_ok,
        $fq_files, $fq_indexes, $indexed
    ) = @_;

    die "$current_version\n" if $version;

    print_usage()
      and die "WARNING: To run successfully, remember to remove '--help'.\n"
      if $help;
    print_usage()
      and die "ERROR: Missing barcode or path to barcode file.\n"
      unless defined $barcode;
    print_usage()
      and die "ERROR: Missing Sample/Experiment ID ('--id').\n"
      unless defined $id;
    print_usage()
      and die "ERROR: Allowed number of mismatched ('--mismatches') must be a positive integer.\n"
      unless $mismatches_ok =~ /^\d+$/;
    print_usage()
      and die "ERROR: Missing path to FASTQ file(s).\n"
      unless scalar @$fq_files > 0;
    print_usage()
      and die "ERROR: Number of FASTQ files must equal number of index files when using '--indexed'.\n"
      if $indexed && scalar @$fq_files != scalar @$fq_indexes;
}

sub print_usage {
    my $prog = basename($0);

    warn <<"EOF";

USAGE
  $prog [options] -b BARCODE IN.FASTQ

DESCRIPTION
  Extracts fastq reads for specified barcode(s) from one or multiple FASTQ files.
  Use wildcards ('*') to match multiple input FASTQ files.

OPTIONS
  -h, --help                 Print this help message
  -v, --version              Print version number
  --id                       Sample or Experiment ID
  -b, --barcode   BARCODE    Specify barcode or file w/ list of barcodes to extract
  -l, --list                 Indicate BARCODE is a list of barcodes in a file
  --indexed                  Samples designated by index sequences
                              Alternate read FQ files and index FQ files
  -m, --mismatches           Minimum number of mismatches allowed in barcode sequence [0]
  -n, --notrim               Split without trimming barcodes
  -st, --stats               Output summary stats only (w/o creating fastq files)
  -o, --outdir    DIR        Output file is saved in the specified directory
                              (or same directory as IN.FASTQ, if --outdir is not used)

NAMING OPTIONS
  --autoprefix               Append FASTQ file name onto output
  --autosuffix               Append barcode onto output
  -p, --prefix    PREFIX     Add custom prefix to output
  -su, --suffix   SUFFIX     Add custom suffix to output

OUTPUT
  An output file in fastq format is written for each barcode to the directory
  containing IN.FASTQ, unless an output directory is specified.
  The default name of the output file is ID.fq. The output names can be
  customized using the Naming Options.

EXAMPLES
  $prog -b GACTG -i Charlotte kitten_DNA.fq
  $prog --barcode barcode.file --list *_DNA.fastq
  $prog --help

EOF
}

sub get_barcodes {
    my ( $list, $barcode, $id ) = @_;

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

    return \%barcode_table;
}

sub validate_barcodes {
    my $barcode_table = shift;

    my @barcodes   = keys %{$barcode_table};
    my $min_length = min map {length} @barcodes;
    my $max_length = max map {length} @barcodes;
    die
      "Unexpected variation in barcode length (min=$min_length, max=$max_length)\n"
      unless $min_length == $max_length;
    my $barcode_length = $max_length;
    map { die "Invalid barcode found: $_\n" unless /^[ACGT]{$barcode_length}$/i }
      @barcodes;

    return $barcode_length;
}

sub parse_filenames {
    my ( $fq_files, $outdir, $barcode ) = @_;

    my ( $filename, $directory, $filesuffix )
        = fileparse( $$fq_files[0], ".f(ast)?q" );
    $filename  = "multi_fq"    if @$fq_files > 1;
    $directory = $outdir . "/" if defined $outdir;
    my $barcode_name = fileparse($barcode);

    return $directory, $filename, $barcode_name;
}

sub open_fq_fhs {
    my ( $barcode_table, $directory, $filename, $barcode_name, $prefix, $suffix, $autoprefix, $autosuffix ) = @_;

    $prefix .= "." unless $prefix eq "";
    $suffix = "." . $suffix unless $suffix eq "";
    $prefix = join ".", $filename, $prefix if $autoprefix;

    my $unmatched_fq_out
        = $directory
        . "unmatched."
        . $prefix . "fq_"
        . $filename . ".bar_"
        . $barcode_name
        . $suffix . ".fq";
    open my $unmatched_fh, ">", $unmatched_fq_out;

    for ( keys %{$barcode_table} ) {
        my $temp_suffix = $suffix;
        $temp_suffix = join ".", $suffix, $_ if $autosuffix;
        my $fq_out
            = $directory
            . $prefix
            . $$barcode_table{$_}->{id}
            . $temp_suffix . ".fq";
        open $$barcode_table{$_}->{fh}, ">", $fq_out;
    }

    return $unmatched_fh;
}

sub split_trim_barcodes {
    my ( $fq_files, $fq_indexes, $indexed, $barcode_table, $barcode_length,
        $unmatched_fh, $notrim, $stats )
        = @_;

    my $total_matched   = 0;
    my $total_unmatched = 0;
    my %barcodes_obs;

    for my $i ( 0 .. $#$fq_files ) {

        my $fq_in     = $$fq_files[$i];
        # my $fq_idx_in = $$fq_files[$i] if $indexed;
        my $fq_idx_in = $$fq_indexes[$i] if $indexed;

        open my $fq_in_fh,     "<", $fq_in;
        open my $fq_idx_in_fh, "<", $fq_idx_in if $indexed;

        while ( my $read_id = <$fq_in_fh> ) {
            my $line_no = $.;
            my $seq     = <$fq_in_fh>;
            my $qual_id = <$fq_in_fh>;
            my $qual    = <$fq_in_fh>;
            validate_fq_read( $read_id, $seq, $qual, $fq_in, $line_no);

            my $cur_barcode = substr $seq, 0, $barcode_length;

            if ($indexed) {
                my $idx_read_id = <$fq_idx_in_fh>;
                my $idx_line_no = $.;
                my $idx_seq     = <$fq_idx_in_fh>;
                my $idx_qual_id = <$fq_idx_in_fh>;
                my $idx_qual    = <$fq_idx_in_fh>;
                validate_fq_read(
                    $idx_read_id, $idx_seq, $idx_qual,
                    $fq_idx_in,   $idx_line_no
                );

                chomp $idx_seq;
                $cur_barcode = $idx_seq;
            }

            my $has_match = exists $$barcode_table{$cur_barcode} ? 1 : 0;

            ( $cur_barcode, $has_match )
                = fuzzy_match( $cur_barcode, $barcode_table, $mismatches_ok )
                if !$has_match && $mismatches_ok;

            $barcodes_obs{$cur_barcode}++;
            if ($has_match) {
                if ( !$indexed && !$notrim ) {
                    $seq  = substr $seq,  $barcode_length + 1;
                    $qual = substr $qual, $barcode_length + 1;
                }
                print { $$barcode_table{$cur_barcode}->{fh} }
                  $read_id . $seq . $qual_id . $qual
                  unless $stats;
                $$barcode_table{$cur_barcode}->{count}++;
                $total_matched++;
            }
            else {
                print $unmatched_fh $read_id . $seq . $qual_id . $qual
                  unless $stats;
                $total_unmatched++;
            }
        }

        close $fq_in_fh;
        close $fq_idx_in_fh if $indexed;
    }

    return $total_matched, $total_unmatched, \%barcodes_obs;
}

sub validate_fq_read {
    my ( $read_id, $seq, $qual, $fq_in, $line_no ) = @_;

    die chomp $read_id
        && "Encountered sequence ID ($read_id) that doesn't start with '\@' on line $line_no of FASTQ file: '$fq_in'...\nInvalid or corrupt FASTQ file?\n"
        unless $read_id =~ /^@/;
    die chomp $read_id
        && "Encountered unequal sequence and quality lengths for read ($read_id) near line $line_no of FASTQ file: '$fq_in'...\nInvalid or corrupt FASTQ file?\n"
        unless length $seq == length $qual;
    die chomp $read_id && chomp $seq && "Encountered sequence ($seq) containing non-nucleotide characters. See read ($read_id) starting on line $line_no of FASTQ file: '$fq_in'...\nInvalid or corrupt FASTQ file?\n"
        unless $seq =~ /^[ACGTN]+$/i;
}

sub fuzzy_match {
    my ( $cur_barcode, $barcode_table, $mismatches_ok ) = @_;

    my %edit_distances;
    push @{ $edit_distances{ distance( $cur_barcode, $_ ) } }, $_
        for keys %{$barcode_table};

    my $best_score       = min keys %edit_distances;
    my $best_score_count = scalar @{ $edit_distances{$best_score} };

    my $fuzzy_barcode;
    my $fuzzy_match;

    if ( $best_score <= $mismatches_ok && $best_score_count == 1 ) {
        ($fuzzy_barcode) = @{ $edit_distances{$best_score} };
        $fuzzy_match = 1;
    }
    else {
        $fuzzy_barcode = $cur_barcode;
        $fuzzy_match   = 0;
    }

    return $fuzzy_barcode, $fuzzy_match;
}

sub close_fq_fhs {
    my ( $barcode_table, $unmatched_fh ) = @_;

    map { close $$barcode_table{$_}->{fh} } keys %{$barcode_table};
    close $unmatched_fh;
}

sub summarize_observed_barcodes {
    my ( $barcode_table, $barcodes_obs, $total_count, $directory, $filename, $barcode_name ) = @_;

    my @sorted_barcodes_obs
        = map {
        [   $_->[0],
            commify( $$barcodes_obs{ $_->[0] } ),
            percent( $$barcodes_obs{ $_->[0] } / ($total_count) ),
            $_->[2]->{id},
        ]
        }
        sort { $b->[1] <=> $a->[1] or $a->[0] cmp $b->[0] }
        map { [ $_, $$barcodes_obs{$_}, $$barcode_table{$_} ] }
        keys %{$barcodes_obs};
    open my $bar_log_fh, ">", $directory . join ".", "log_barcodes_observed",
        "fq_" . $filename, "bar_" . $barcode_name;
    my $tbl_observed = Text::Table->new(
        "barcode",
        "count\n&right",
        "percent\n&right",
        "id\n&right",
    );
    $tbl_observed->load(@sorted_barcodes_obs);
    print $bar_log_fh $tbl_observed;
    close $bar_log_fh;
}

sub summarize_counts {
    my ( $barcode_table, $fq_files, $total_count, $total_matched,
        $total_unmatched, $directory, $filename, $barcode_name )
        = @_;

    my @barcode_counts =
      sort { $a->[0] cmp $b->[0] }
      map { [ $$barcode_table{$_}->{id}, $_, commify( $$barcode_table{$_}->{count} ), percent( $$barcode_table{$_}->{count} / $total_count ) ] }
      keys %{$barcode_table};

    open my $count_log_fh, ">", $directory . join ".", "log_barcode_counts", "fq_" . $filename, "bar_" . $barcode_name;
    say $count_log_fh "Barcode splitting summary for:";
    map { say $count_log_fh "  " . $_ } @$fq_files;
    my $max_fq_length = max map { length } @$fq_files;
    say $count_log_fh "-" x ( $max_fq_length + 2 );

    my $tbl_matched = Text::Table->new(
        "",
        "\n&right",
        "\n&right",
    );
    $tbl_matched->load(
        [
            "matched",
            commify($total_matched),
            percent( $total_matched / ( $total_matched + $total_unmatched ) ),
        ],
        [
            "unmatched",
            commify($total_unmatched),
            percent( $total_unmatched / ( $total_matched + $total_unmatched ) ),
        ],
    );

    print $count_log_fh $tbl_matched;
    say $count_log_fh "-" x ( $max_fq_length + 2 );

    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data( map{ $$barcode_table{$_}->{count} } keys %{$barcode_table} );
    my $tbl_stats = Text::Table->new(
        "",
        "\n&num",
        "\n&right",
    );
    $tbl_stats->load(
        [ "barcodes", commify( $stat->count ) ],
        [ "min",      commify( $stat->min ),           percent( $stat->min / $total_count ) ],
        [ "max",      commify( $stat->max ),           percent( $stat->max / $total_count ) ],
        [ "mean",     commify( round( $stat->mean ) ), percent( round( $stat->mean ) / $total_count ) ],
        [ "median",   commify( $stat->median ),        percent( $stat->median / $total_count ) ],
    );
    print $count_log_fh $tbl_stats;
    say $count_log_fh "-" x ( $max_fq_length + 2 );

    my $tbl_counts = Text::Table->new(
        "id",
        "barcode",
        "count\n&right",
        "percent\n&right",
    );
    $tbl_counts->load(@barcode_counts);
    print $count_log_fh $tbl_counts;
    close $count_log_fh;
}

sub plot_summary {
    my ( $barcodes_obs, $barcode_table, $directory, $id ) = @_;

    my ( $barcode_data, $count_data, $match_data )
        = get_vectors( $barcodes_obs, $barcode_table );

    my $R = Statistics::R->new();

    $R->run(qq`setwd("$directory")`);
    $R->run(qq`log <- data.frame(barcode = $barcode_data,
                                 count   = $count_data,
                                 matched = $match_data)`);

    my $plot_fxn = plot_function();

    $R->run($plot_fxn);
    $R->run(qq`barcode_plot(log, "$id")`);
}

sub get_vectors {
    my ( $barcodes_obs, $barcode_table ) = @_;

    my @barcodes;
    my @counts;
    my @matches;
    for my $barcode ( keys %{$barcodes_obs} ) {
        push @barcodes, $barcode;
        push @counts, $$barcodes_obs{$barcode};
        my $match = exists $$barcode_table{$barcode} ? "matched" : "unmatched";
        push @matches, $match;
    }

    my $barcode_data = join ',', map {qq!"$_"!} @barcodes;
    my $count_data   = join ',', @counts;
    my $match_data   = join ',', map {qq!"$_"!} @matches;

    $_ = "c($_)" for $barcode_data, $count_data, $match_data;

    return $barcode_data, $count_data, $match_data;
}

sub plot_function {
    return <<EOF;    # Adapted from barcode_plot.R
        barcode_plot <- function(log.df, id) {
          library("ggplot2")

          log.plot <-
            ggplot(
              data = log.df,
              aes(x = factor(matched), y = count / 1000000)) +
            geom_boxplot(outlier.size = 0) +
            geom_jitter() +
            ggtitle(label = id) +
            xlab(label = "") +
            scale_y_continuous(name = "Count\n(million reads)")

          ggsave(
            filename = paste(sep = "", id, ".barcodes.png"),
            plot     = log.plot,
            width    = 4,
            height   = 5
          )
        }
EOF
}
