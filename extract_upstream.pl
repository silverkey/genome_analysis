#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

use Bio::DB::SeqFeature::Store;
use Bio::Seq;
use Bio::SeqIO;

my $USAGE = "\nperl $0 [database name] [start of the fragment (number of bases before TSS)] [end of the fragment (number of bases before TSS)] [gene|gene:protein_coding]\n\n";

die $USAGE unless $ARGV[0];
die $USAGE unless $ARGV[1];
die $USAGE unless $ARGV[2];
die $USAGE unless $ARGV[3];
die $USAGE unless ($ARGV[3] eq 'gene' || $ARGV[3] eq 'gene:protein_coding');

my $dbname = $ARGV[0];
my $START = $ARGV[1];
my $END = $ARGV[2];
my $TYPE = $ARGV[3];

my $MINIMUM_LENGTH = $START - $END - 10;
my $MAXIMUM_LENGTH = $START - $END + 10;

die "\n Start must be bigger than end!\n\n" unless $START > $END;
die "\n Start and end must be positive numbers!\n\n" unless $START > 0 && $END > 0;

my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                         -user => 'mysql_dev',
                                         -pass => 'riiGbs',
                                         -dsn => "dbi:mysql:$dbname");

my $tab_out = "$dbname\_upstream_$START\_$END\.tab";
my $fa_out = "$dbname\_upstream_$START\_$END\.fa";

open(OUT,">$tab_out");
my $seqio = Bio::SeqIO->new(-file => ">$fa_out",
                            -format => 'fasta');

my @genes = $db->get_features_by_type($TYPE);

foreach my $gene(@genes) {

  my $seq_id = $gene->seq_id;
  my @attributes = $gene->attributes('load_id');
  my $load_id = $attributes[0];
  my $id = $gene->id || 'NA';
  my $source = $gene->source;
  my $strand = $gene->strand;
  my $start = $gene->start;
  my $end = $gene->end;

  my $pstart;
  my $pend;
  if($strand == 1) {
    $pstart = $start-$START;
    $pend = $start-$END;
  }
  elsif($strand == -1) {
    $pstart = $end+$END;
    $pend = $end+$START;
  }

  my $desc = "$id $seq_id\:$pstart\-$pend $strand";
  
  my $string = $db->fetch_sequence(-seq_id => $seq_id,
                                   -start => $pstart,
                                   -end => $pend);

  my $seq = Bio::Seq->new(-id => $load_id,
                          -desc => $desc,
                          -seq => $string);

  if($strand == -1) {
    my $revcom = $seq->revcom;
    $seq->seq($revcom->seq);
  }

  my $pseq = $seq->seq;

  if($seq->length >= $MINIMUM_LENGTH && $seq->length <= $MAXIMUM_LENGTH) {
    print OUT "$load_id\t$id\t$seq_id\t$start\t$end\t$strand\t$pstart\t$pend\t$pseq\n" if $seq->length > 10;
    $seqio->write_seq($seq) if $seq->length > 10;
  }
}
