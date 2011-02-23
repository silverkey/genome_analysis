#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use Data::Dumper;

use Bio::DB::SeqFeature::Store;
use Bio::Seq;
use Bio::SeqIO;

use lib '/home/remo/src/svn/UTILS/';
use Range::Utils;

my $usage = "\n\tperl $0 [db] [type (e.g.: exon or exon:protein_coding)]\n\n";

die $usage if scalar(@ARGV) != 2;
die $usage unless ($ARGV[1] eq 'exon' || $ARGV[1] eq 'exon:protein_coding');

my $ID = 1;
my $type = $ARGV[1];

# FILES AND DATABASES:
my $dbname = $ARGV[0];
my $faname = "$dbname\_exonic.fa";

my $seqio = Bio::SeqIO->new(-file => ">$faname",
                            -format => 'fasta');

# DATABASES CONNCTIONS:
# TEST DATABASE TO PLAY WITH RANGES
my $dbtest = DBI->connect('dbi:mysql:test','mysql_dev','riiGbs',{PrintWarn=>1,PrintError=>1,RaiseError=>1}) or die $DBI::errstr;
# SeqFeature DATABASE OF THE GENOME ANNOTATIONS
my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                         -user => 'mysql_dev',
                                         -pass => 'riiGbs',
                                         -dsn => "dbi:mysql:$dbname");

my @seq_ids = $db->seq_ids;

foreach my $seq_id(@seq_ids) {

  my $seq = $db->fetch_sequence(-seq_id => $seq_id,
                                -bioseq => 1);

  my $gene_ranges = get_gene_ranges($db,$seq);

  if($gene_ranges) {
    my $disconnected = Range::Utils->disconnected_ranges($dbtest,$gene_ranges,1);
    extract_sequences($disconnected,$seq);
  }

  else {
    # IF WE HAVE NO GENES ACTUALLY WE DO NOT CONSIDER THE SCAFFOLD
    # toplevel_with_no_genes($seq);
  }
}

sub get_gene_ranges {
  my $db = shift;
  my $seq = shift;
  print "Working on ".$seq->id."\n";

  my $aref;

  my @genes = $db->features(-seq_id => $seq->id,
                            -type => $type);

  foreach my $gene(@genes) {
    my $range = Bio::Range->new(-start => $gene->start,
                                -end => $gene->end);
    push(@$aref,$range);
  }
  return $aref;
}

sub toplevel_with_no_genes {
  my $seq = shift;
  my $start = 1;
  my $end = $seq->length;

  annotate_fragment($seq,$start,$end) if $end - $start > 10;
}

sub extract_sequences {
  my $ranges = shift;
  my $seq = shift;
  foreach my $r(@$ranges) {
    my $start = $r->start;
    my $end = $r->end;
    annotate_fragment($seq,$start,$end) if $end - $start > 10;
  }
}

sub annotate_fragment {
  my $seq = shift;
  my $start = shift;
  my $end = shift;

  my $last = Bio::Seq->new(-id => $ID,
                           -desc => $seq->id.':'.$start.'-'.$end,
                           -seq => $seq->subseq($start,$end));
  $seqio->write_seq($last);

  $ID++;
}
