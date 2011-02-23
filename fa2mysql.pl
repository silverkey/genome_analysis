#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;
use Bio::SeqFeature::Generic;

my $USAGE = "\nperl $0 [fasta file] [dbname] [create 'Y' or 'N']\n\n";

die $USAGE unless -e $ARGV[0];
die $USAGE unless $ARGV[1];
die $USAGE unless $ARGV[2];

my $fa_file = $ARGV[0];
my $dbname = $ARGV[1];
my $create = 0;
$create = 1 if $ARGV[2] eq 'Y';

my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                         -user => 'mysql_dev',
                                         -pass => 'riiGbs',
                                         -dsn => "dbi:mysql:$dbname",
                                         -create => $create);

my $seqio = Bio::SeqIO->new(-format => 'fasta',
                            -file => $fa_file);

my @features;

while(my $seq = $seqio->next_seq) {
  # --------------------------------------------------------------------------------
  # NOT GENERAL LINES FOR DIATOMS DATA DOWNLOADED FROM THE ENSEMBL!!!
  # --------------------------------------------------------------------------------
  my $id = $seq->id;
  if($id eq 'EG:') {
    # This is because same ID from the ensembl fasta look like the next row line:
    # EG: bd_34x588 dna_rm:supercontig supercontig:Thaps3_bd:bd_34x588:1:2282:1
    # Therefore we will need to use the name after the first spacein order to get the ID
    $id = $seq->description;
    $id =~ s/([\d\w]+) .+/$1/;
  }
  $id =~ s/EG\://;
  print '-'.$id."-\n";
  $db->insert_sequence($id,$seq->seq);

  my $feature = Bio::SeqFeature::Generic->new(-seq_id => $id,
                                              -start => 1,
                                              -end => $seq->length,
                                              -primary => 'toplevel',
                                              -source => 'JGI',
                                              -strand => 1,
                                              -attributes => { ID => $id } );
  push(@features,$feature);
}

$db->store(@features) or die "\nPROBLEM!!!\n";
