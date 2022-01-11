#!/usr/bin/perl
# script imported from 
# https://github.com/lindberg-m/count_kmers

use strict;
use warnings;

my $usage = << "EOF";
usage:
  $0 [-h,-v,-p,-i,-u,-l,-b bedfile, -m bedfile] [SIZE] < fasta_file.fa
  Optional Arguments
    -h, --help        Show this message and exit.
    -v, --verbose     Increase terminal output.
    -p, --pyrimidine  When kmers are of odd length, use pyrimidine based
                      contexts with the middle nucleotide as a reference point
    -i, --ignore-amb  Ignore ambiguous contexts [default=true], setting this
                      param turn this feature off.
    -u, --upper       Turn all sequences uppercase before counting.
    -l, --ignore-low  Ignore contexts with lower case characters in them.
    -b, --bed  FILE   Bedfile that indicate what regions to perform counting on.
    -m, --mask FILE   Bedfile to use for masking fasta. Positions overlapping
                      mask regions will be transformed to 'n' characters.
    -s  --size SIZE   Size of kmers to calculate bgfreq on. Defaults to 3

This program reads a fasta formatted file from standard input and
calculate k-mer frequencies of specified length. Optionally,
subset or mask sequences based on bedfile(s) first.
EOF

my %PARAMS = (                  
  SIZE       => 3,  # K-mer size
  VERBOSE    => 0,  # Terminal Verbosity
  PYRIMIDINE => 0,  # Report pyrimidine based kmers only
  IGNORE_AMB => 1,  # Remove ambiguous nucleotides
  MK_UPPER   => 0,  # Make sequences uppercase before any processing
  IGNORE_LOW => 0,  # Remove lower case nucleotides
  BEDFILE    => '', # Path to bed-file for regions to do analysis on
  MASK       => '', # Path to bed-file to use for masking the genome
);
  

sub main {
  parse_args(\%PARAMS);
  my $kmer_counts = {};
  my $regions     = $PARAMS{BEDFILE} ? parse_bed($PARAMS{BEDFILE}) : {};
  my $mask        = $PARAMS{MASK}    ? parse_bed($PARAMS{MASK})    : {};

  # Parse fasta from STDIN and perform subsetting and counting
  my $last_chrom = '';
  my $seq        = '';
  while (<STDIN>) {
    chomp;
    if ( /^>/ ) {
      s/^>//; print STDERR "$_\n" if ($PARAMS{VERBOSE});
      count_kmers(\$seq, $mask->{$last_chrom}, $regions->{$last_chrom}, $kmer_counts) if ($seq);
      $last_chrom = $_;
      $seq = '';
    } else {
      $seq .= $_;
    }
  }

  # Last sequence isn't handeled in the while-loop, hence this:
  count_kmers(\$seq, $mask->{$last_chrom}, $regions->{$last_chrom}, $kmer_counts) if ($seq);

  # Combine pyrimidine based kmers with its reverse complement counterparts if necessary
  my $kmer_count_results = {};
  if ($PARAMS{PYRIMIDINE}) {
    my $mid_point = ($PARAMS{SIZE} - 1) / 2;
    for my $k (keys %{$kmer_counts}) {
      my $k2 = (substr($k, $mid_point, 1) =~ /cCtT/) ? $k : revcomp($k);
      $kmer_count_results->{$k2} += $kmer_counts->{$k};
    }
  } else {
    $kmer_count_results = $kmer_counts;
  }

  # Then print results
  for my $k (sort keys %{$kmer_count_results}) {
    print "$k\t$kmer_count_results->{$k}\n";
  }
}

sub count_kmers {
  my $seq = shift;     # Scalar ref to dna sequence
  my $mask = shift;    # Arref      to regions that should be masked
  my $regions = shift; # Arref      to regions to performa analysis on
  my $counts = shift;  # Hashref    to kmer-counts

  mask_sequence($seq, $mask) if ($PARAMS{MASK});
  $$seq = uc $$seq           if ($PARAMS{MK_UPPER});

  if ($PARAMS{BEDFILE}) {
    foreach my $region (@{$regions}) {

      my $start  = $region->[0];
      my $end    = $region->[1];
      my $subseq = substr($$seq, $start, $end - $start);

      my $seqparts = remove_ambig_and_softmasked($subseq);
      map { update_counts($_, $counts) } @{$seqparts};
    }
  } else {
    my $seqparts = remove_ambig_and_softmasked($$seq);
    map { update_counts($_, $counts) } @{$seqparts};
  }
}

sub remove_ambig_and_softmasked {
  my $seq          = shift; # Sequence

  my (@temp, @res);
  if ($PARAMS{IGNORE_LOW}) {
    foreach (split /[a-z]+/, $seq) {
      push @temp, $_;
    }
  } else {
    push @temp, $seq;
  }

  if ($PARAMS{IGNORE_AMB}) {
    foreach my $sp (@temp) {
      foreach (split /[nxywrNXYWR]+/, $sp) {
        push @res, $_;
      }
    }
  }
  return \@res;
}

sub mask_sequence {
  my $seqr = shift; # Scalar Ref to sequence
  my $mask = shift; # Arref to (start, end) pairs

  for my $r (@$mask) {
    my $seqlen = $r->[1] - $r->[0];
    substr($$seqr, $r->[0], $seqlen) = 'n' x $seqlen;
  }
}

sub update_counts {
  my $seq       = shift; # Scalar,  DNA sequences
  my $counts    = shift; # Hashref, kmer counts

  my ($ctx, $start, $end);
  my $seqlen = length($seq);
  for ($start=0,$end=$PARAMS{SIZE}; $end <= $seqlen; $end++,$start++) {
    $ctx = substr($seq, $start, $PARAMS{SIZE});
    $counts->{$ctx}++;
  }
}

sub parse_bed {
  my $bedfile = shift;
  my %regions;
  open BED, '<', $bedfile or die "Cannot open bedfile: $bedfile\n";
  while (<BED>) {
    chomp;
    my ($chrom, $start, $stop) = split;
    push @{$regions{$chrom}}, [ $start, $stop ];
  }
  return \%regions;
}

sub parse_args {
  my $pars = shift;
  my $j = 0;
  for (my $i = 0; $i<@ARGV; $i++) {
    if ($ARGV[$i] =~ /^-/) {
      if ($ARGV[$i] eq '-h' || $ARGV[$i] eq '--help') {
        die $usage;
      } elsif ($ARGV[$i] eq '-v' || $ARGV[$i] eq '--verbose') {
        $pars->{VERBOSE} = 1;
      } elsif ($ARGV[$i] eq '-p' || $ARGV[$i] eq '--pyrimidine') {
        $pars->{PYRIMIDINE} = 1;
      } elsif ($ARGV[$i] eq '-i' || $ARGV[$i] eq '--ignore-amb'){
        $pars->{IGNORE_AMB} = 0;
      } elsif ($ARGV[$i] eq '-u' || $ARGV[$i] eq '--upper'){
        $pars->{MK_UPPER} = 1;
      } elsif ($ARGV[$i] eq '-l' || $ARGV[$i] eq '--ignore-low'){
        $pars->{IGNORE_LOW} = 1;
      } elsif ($ARGV[$i] eq '-b' || $ARGV[$i] eq '--bed') {
        $pars->{BEDFILE} = $ARGV[++$i];
      } elsif ($ARGV[$i] eq '-m' || $ARGV[$i] eq '--mask') {
        $pars->{MASK} = $ARGV[++$i]
      } elsif ($ARGV[$i] eq '-s' || $ARGV[$i] eq '--size') {
        $pars->{SIZE} = $ARGV[++$i]      
      } else {
        die "Unrecongnized argument: $ARGV[$i]\n";
      }
    } else {
      # $j++;
      # if ($i == scalar @ARGV) {
        # $pars->{SIZE} = $ARGV[$i];
      #do nothing
      ;
    }
  }

  die "Need at most 1 positional argument but got $j\n" . $usage if ($j > 1);
  if ($pars->{SIZE} % 2 == 0 && $pars->{PYRIMIDINE}) {
    die "Cannot use pyrimidine based calculation on an even-sized k-mers\n";
  }
}

sub revcomp {
  my $sequence = shift;
  $sequence =~ tr/actgACTG/tgacTGAC/;
  return reverse $sequence;
}

main()