#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $upfile = 'FC_upstream_100_1.fa';
my $downfile = 'FC_downstream_1_100.fa',
my $word = 'CGACGA';

my $up = Bio::SeqIO->new(-file => $upfile,
                         -format => 'fasta');

my $down = Bio::SeqIO->new(-file => $downfile,
                           -format => 'fasta');

my $href;
my $totup;
my $totdown;
my $upn;
my $downn;
my $common;

while(my $seq = $up->next_seq) {
  my $string = uc($seq->seq);
  if($string =~ /$word/) {
    $href->{$seq->id}->{up} ++;
    $upn ++;
  }
  $totup ++;
}

while(my $seq = $down->next_seq) {
  my $string = uc($seq->seq);
  if($string =~ /$word/) {
    $href->{$seq->id}->{down} ++;
    $downn ++;
  }
  $totdown ++;
}

foreach my $id(keys %$href) {
  next unless exists $href->{$id}->{up};
  next unless exists $href->{$id}->{down};
  $common ++;
}

print "Total up sequences containing $word: $totup\n".
      "Total down sequences containg $word: $totdown\n".
      "Upstream containing $word: $upn\n".
      "Downstream containing $word: $downn\n".
      "Upstream and downstream: $common\n".
      "Frequency in up: ".($upn/$totup)."\n".
      "Frequency in down: ".($downn/$totdown)."\n".
      "Expected upstrea and downstream: ".($upn/$totup)*($downn/$totdown)*$totdown."\n";
