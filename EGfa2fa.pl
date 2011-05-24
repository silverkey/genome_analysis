#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $USAGE = "\nperl $0 [fasta file]\n\n";

die $USAGE unless -e $ARGV[0];

my $EG_fa_file = $ARGV[0];

my $fa_file = $ARGV[0].'_mod';

my $seqin = Bio::SeqIO->new(-format => 'fasta',
                            -file => $EG_fa_file);

my $seqout = Bio::SeqIO->new(-format => 'fasta',
                             -file => ">$fa_file");

while(my $seq = $seqin->next_seq) {
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

  $seq->id($id);
  $seqout->write_seq($seq);
}

