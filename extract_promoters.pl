#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

use Bio::DB::SeqFeature::Store;
use Bio::Seq;
use Bio::SeqIO;

my $USAGE = "\nperl $0 [database name] [pretss length] [posttss length] [type gene|gene:protein_coding]\n\n";

die $USAGE unless $ARGV[0];
die $USAGE unless $ARGV[1];
die $USAGE unless $ARGV[2];
die $USAGE unless $ARGV[3] eq 'gene' || $ARGV[3] eq 'gene:protein_coding';

my $dbname = $ARGV[0];
my $pre_tss = $ARGV[1];
my $post_tss = $ARGV[2];
my $type = $ARGV[3];

my $MINIMUM_LENGTH = $post_tss + $pre_tss - 10;
my $MAXIMUM_LENGTH = $post_tss + $pre_tss + 10;

my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                         -user => 'mysql_dev',
                                         -pass => 'riiGbs',
                                         -dsn => "dbi:mysql:$dbname");

my $tab_out = "$dbname\_promoter_$pre_tss\_$post_tss\.tab";

open(OUT,">$tab_out");
my $fa_out = "$dbname\_promoter_$pre_tss\_$post_tss\.fa";
my $seqio = Bio::SeqIO->new(-file => ">$fa_out",
                            -format => 'fasta');

my @genes = $db->get_features_by_type($type);

foreach my $gene(@genes) {

  my $seq_id = $gene->seq_id;
  my @attributes = $gene->attributes('load_id');
  my $load_id = $attributes[0];
  my $id = $gene->id || 'NA';
  my $source = $gene->source;
  my $strand = $gene->strand;
  my $start = $gene->start;
  my $end = $gene->end;

#  next unless $source eq 'protein_coding';

  my $pstart;
  my $pend;
  if($strand == 1) {
    $pstart = $start-$pre_tss;
    $pend = $start+$post_tss;
  }
  elsif($strand == -1) {
    $pstart = $end-$post_tss;
    $pend = $end+$pre_tss;
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
