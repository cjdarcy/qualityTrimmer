#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

# Author: Christopher Darcy
# qualityTrimmer.pl: Combines left and right FASTQ-formatted
# reads into an interleaved FASTQ file with a given
# quality score threshold.

# Modules
use Bio::Seq;
use Bio::Seq::Quality;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

# Globals
my $left	= '';
my $right	= '';
my $interleaved = '';
my $qual	= 0;
my $usage	= "\n$0 [options] \n
Options:
	-left		Left reads
	-right		Right reads
	-qual		Quality score minimum
	-interleaved	Filename for interleaved output
	-help		Show this message
\n";

# Check the flags
GetOptions(
	'left=s'	=> \$left,
	'right=s'	=> \$right,
	'interleaved=s' => \$interleaved,
	'qual=i'	=> \$qual,
	'help'		=> sub { pod2usage($usage); },
) or pod2usage($usage);

# Error checking for options
unless (-e $left and -e $right and $qual and $interleaved ) {
	unless (-e $left) {
		print "Specify file for left reads\n";
	}
	unless (-e $right) {
		print "Specify file for right reads\n";
	}
	unless ($interleaved) {
		print "Specify file for interleaved output\n";
	}
	unless ($qual) {
		print "Specify quality score cutoff\n", $usage;
	}
	die "Missing required options\n";
}

# Create output file
my $seq_out = Bio::SeqIO->new(
        -file => ">$interleaved",
        -format => 'fastq'
        ) or die $!;

# Open paired end reads (R1 and R2)
my $R1 = Bio::SeqIO->new(
	-file => $left,
	-format => 'fastq'
	) or die $!;

my $R2 = Bio::SeqIO->new(
	-file => $right,
	-format => 'fastq'
	) or die $!;

# Loop through input files
# Return longest subsequences above a given quality score cutoff
# Write to interleaved output file
while ( my $seq_obj_1 = $R1->next_seq() ) {
	my $seq_obj_2 = $R2->next_seq();
	my $leftTrimmed = $seq_obj_1->get_clear_range($qual);
	my $rightTrimmed = $seq_obj_2->get_clear_range($qual);
	$leftTrimmed->desc($seq_obj_1->desc());
	$rightTrimmed->desc($seq_obj_2->desc());
	$seq_out->write_seq($leftTrimmed);
	$seq_out->write_seq($rightTrimmed);
}
