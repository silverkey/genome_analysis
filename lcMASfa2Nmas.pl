#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

my $USAGE = "\nperl $0 [fasta file]\n\n";

die $USAGE unless -e $ARGV[0];

my $fa_file = $ARGV[0];

my $ma_fa_file = $ARGV[0].'_N_mask';

my $seqin = Bio::SeqIO->new(-format => 'fasta',
                            -file => $fa_file);

my $seqout = Bio::SeqIO->new(-format => 'fasta',
                             -file => ">$ma_fa_file");

while(my $seq = $seqin->next_seq) {
  my $seq = $seq->seq;
  $seq =~ s/[atcg]/N/g;
  $seq->seq($seq);
  $seqout->write_seq($seq);
}
