# Barcode Splitter, Trimmer, and Stat Generator

Contents of README:

- [Usage Summary][usage_summary]
- [Detailed Usage Example][detailed_example]
    - [Barcode splitting/trimming][split_trim]
    - [Plotting barcode splitting summary][plot]

## Usage Summary

The following can be accessed by running `./barcode_split_trim.pl --help`:

    USAGE
      barcode_split_trim.pl [options] -b BARCODE IN.FASTQ

    DESCRIPTION
      Extracts fastq reads for specified barcode(s) from one or multiple FASTQ files.
      Use wildcards ('*') to match multiple input FASTQ files.

    OPTIONS
      -h, --help                 Print this help message
      -v, --version              Print version number
      -b, --barcode   BARCODE    Specify barcode or list of barcodes to extract
      -i, --id        SAMPLE_ID  Sample ID (not needed if using list of barcodes)
      -l, --list                 Indicates --barcode is a list of barcodes in a file
      -n, --notrim               Split without trimming barcodes off
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
      The default name of the output file is SAMPLE_ID.fq. The output names can be
      customized using the Naming Options.

    EXAMPLES
      barcode_split_trim.pl -b GACTG -i Charlotte kitten_DNA.fq
      barcode_split_trim.pl --barcode barcode.file --list *_DNA.fastq
      barcode_split_trim.pl --help

## Detailed Usage Example

A pair of sample FASTQ and barcode files is provided in the [`sample_files`][sample_files_folder] folder. They have been used to generate the output FASTQ, log, and summary plot files ([`sample_files/output`][output_folder]). In this example, I show what these files look like and how to generate them.

### Barcode splitting/trimming

The sample FASTQ file (`sample_files/sequences.fq`) contains 100,000 sequence reads from a pool of 14 samples. The corresponding [barcodes][barcodes] (`sample_files/barcode.list`) are:

    TACGC   marmot1
    ATCGT   marmot2
    ATTCC   marmot3
    CCAGC   marmot4
    GATAC   aardvark1
    GGATG   aardvark2
    TCGAT   tarsier1
    AGCGC   tarsier2
    CCAAT   tarsier3
    CGCTG   tarsier4
    CTAGT   puffin1
    TAGAG   puffin2
    TAGTC   puffin3
    GGTCA   puffin4

To split `sample_files/sequences.fq` with the barcodes in `sample_files/barcode.list`, we can run:

```bash
./barcode_split_trim.pl \
  --barcode sample_files/barcode.list \
  --list \
  --outdir sample_files/output \
  sample_files/sequences.fq
```

<!--
    # The above snippet is Github's syntax-highlighted version of this:
    ./barcode_split_trim.pl \
      --barcode sample_files/barcode.list \
      --list \
      --outdir sample_files/output \
      sample_files/sequences.fq
 -->

This results in a FASTQ file for each barcode (barcodes are trimmed), a single FASTQ file containing all unmatched barcodes (barcodes are left in tact), and two log files:

    # FASTQ file for each barcode
    aardvark1.fq
    aardvark2.fq
    marmot1.fq
    marmot2.fq
    marmot3.fq
    marmot4.fq
    puffin1.fq
    puffin2.fq
    puffin3.fq
    puffin4.fq
    tarsier1.fq
    tarsier2.fq
    tarsier3.fq
    tarsier4.fq

    # FASTQ file containing all unmatched barcodes
    unmatched.fq_sequences.bar_barcode.list.fq

    # log files
    log_barcode_counts.fq_sequences.bar_barcode.list
    log_barcodes_observed.fq_sequences.bar_barcode.list

The [first log file][LOG1] is the barcode splitting summary (`sample_files/output/log_barcode_counts.fq_sequences.bar_barcode.list`) and returns some basic stats with an emphasis on expected barcodes:

    Barcode splitting summary for:
      sample_files/sequences.fq
    ---------------------------
    matched   24,459 97.8%
    unmatched    541  2.2%
    ---------------------------
    barcodes    14
    min      1,260   5.0%
    max      2,200   8.8%
    mean     1,747   7.0%
    median   1,735.5 6.9%
    ---------------------------
    id        barcode count percent
    aardvark1 GATAC   1,595 6.4%
    aardvark2 GGATG   1,839 7.4%
    marmot1   TACGC   1,758 7.0%
    marmot2   ATCGT   1,504 6.0%
    marmot3   ATTCC   1,614 6.5%
    marmot4   CCAGC   1,468 5.9%
    puffin1   CTAGT   1,476 5.9%
    puffin2   TAGAG   1,260 5.0%
    puffin3   TAGTC   2,102 8.4%
    puffin4   GGTCA   1,964 7.9%
    tarsier1  TCGAT   2,010 8.0%
    tarsier2  AGCGC   2,200 8.8%
    tarsier3  CCAAT   1,713 6.9%
    tarsier4  CGCTG   1,956 7.8%


The [second log file][LOG2] returns counts and percentages for all observed barcodes (both expected and unexpected). Below are the first 20 (of 215) observed barcodes from this log (`sample_files/output/log_barcodes_observed.fq_sequences.bar_barcode.list`):

    barcode count percent id
    AGCGC   2,200 8.8%     tarsier2
    TAGTC   2,102 8.4%      puffin3
    TCGAT   2,010 8.0%     tarsier1
    GGTCA   1,964 7.9%      puffin4
    CGCTG   1,956 7.8%     tarsier4
    GGATG   1,839 7.4%    aardvark2
    TACGC   1,758 7.0%      marmot1
    CCAAT   1,713 6.9%     tarsier3
    ATTCC   1,614 6.5%      marmot3
    GATAC   1,595 6.4%    aardvark1
    ATCGT   1,504 6.0%      marmot2
    CTAGT   1,476 5.9%      puffin1
    CCAGC   1,468 5.9%      marmot4
    TAGAG   1,260 5.0%      puffin2
    GGGCA      32 0.1%
    GGTCC      27 0.1%
    GGATT      15 0.1%
    NAGAG      11 0.0%
    TCTAT      11 0.0%
    NAGTC      11 0.0%


### Plotting barcode splitting summary

The logs are useful, but if there are numerous barcodes and/or experiments being analyzed at once, it can be difficult to easily detect irregularities or problematic barcodes. To solve this problem, we can make a barcode frequency plot using R.

First we extract the [relevant info][LOG2_converted] from the log of observed barcodes:

    ./log_converter.pl sample_files/output/log_barcodes_observed.fq_sequences.bar_barcode.list

Then in R, we define the `barcode_plot` function (found in [`barcode_plot.R`][barcode_plot]) and run the following code:

```R
setwd("sample_files/output/")
expt.name <- "demo"
file.name <- "log_barcodes_observed.fq_sequences.bar_barcode.list.tsv"
barcode_plot(file.name, expt.name)
```

<!--
    # The above snippet is Github's syntax-highlighted version of this:
    setwd("sample_files/output/")
    expt.name <- "demo"
    file.name <- "log_barcodes_observed.fq_sequences.bar_barcode.list.tsv"
    barcode_plot(file.name, expt.name)
 -->

For this plot (saved to `sample_files/output/demo.barcodes.png`), barcodes are split into two groups, those that match an expected barcode and those that are unmatched. Boxplots are then generated using the observed barcode frequencies (which are jitter-plotted individually on top of the boxplot).

<img src="https://raw.github.com/mfcovington/auto_barcode/master/sample_files/output/demo.barcodes.png" height="500" />

<!-- LINKS -->

[usage_summary]: https://github.com/mfcovington/auto_barcode#usage-summary
[detailed_example]: https://github.com/mfcovington/auto_barcode/#detailed-usage-example
[split_trim]: https://github.com/mfcovington/auto_barcode/#barcode-splittingtrimming
[plot]: https://github.com/mfcovington/auto_barcode/#plotting-barcode-splitting-summary

[sample_files_folder]: https://github.com/mfcovington/auto_barcode/tree/master/sample_files
[barcodes]: https://raw.github.com/mfcovington/auto_barcode/master/sample_files/barcode.list

[output_folder]: https://github.com/mfcovington/auto_barcode/tree/master/sample_files/output
[LOG1]: https://raw.github.com/mfcovington/auto_barcode/master/sample_files/output/log_barcode_counts.fq_sequences.bar_barcode.list
[LOG2]: https://raw.github.com/mfcovington/auto_barcode/master/sample_files/output/log_barcodes_observed.fq_sequences.bar_barcode.list
[LOG2_converted]: https://raw.github.com/mfcovington/auto_barcode/master/sample_files/output/log_barcodes_observed.fq_sequences.bar_barcode.list.tsv

[barcode_plot]: https://github.com/mfcovington/auto_barcode/blob/master/barcode_plot.R

