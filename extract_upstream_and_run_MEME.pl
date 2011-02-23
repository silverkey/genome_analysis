#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

use Bio::DB::SeqFeature::Store;
use Bio::Seq;
use Bio::SeqIO;

my $dbname = $ARGV[0]; #  'Phaeodactylum_tricornutum.Phatr2.dna_rm.toplevel.sqlite';

my $START = $ARGV[1]; # 400;
my $END = $ARGV[2]; # 400;
my $BP = $ARGV[3];

my $USAGE = "\nperl $0 [database name] [start of the fragment (number of bases before TSS)] [end of the fragment (number of bases before TSS)] [bp length of motifs for meme]\n\n";

die $USAGE unless $ARGV[0];
die $USAGE unless $ARGV[1];
die $USAGE unless $ARGV[2];
die $USAGE unless $ARGV[3];
die "\n Start must be bigger than end!\n\n" unless $START > $END;
die "\n Start and end must be positive numbers!\n\n" unless $START > 0 && $END > 0;

my $MINIMUM_LENGTH = $START - $END - 10;
my $MAXIMUM_LENGTH = $START - $END + 10;

my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                         -user => 'mysql_dev',
                                         -pass => 'riiGbs',
                                         -dsn => "dbi:mysql:$dbname");

my $tab_out = "$dbname\_upstream_$START\_$END\.tab";
my $fa_out = "$dbname\_upstream_$START\_$END\.fa";

my $dir = "MEME_$fa_out\_$bp";
mkdir($dir);
system("cd $dir");

open(OUT,">$tab_out");
my $seqio = Bio::SeqIO->new(-file => ">$fa_out",
                            -format => 'fasta');

my @genes = $db->get_features_by_type('gene');

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

$seqio->close;
run_meme($fa_out,$BP);

sub run_meme {
  my $file = shift;
  my $bp = shift;
  sytem("nohup meme $file -dna -evt 0.0001 -w $bp -minw $bp -maxw $bp -maxsize 30000000 &");
}
