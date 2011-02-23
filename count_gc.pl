#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;

my $USAGE = "\nUSAGE: perl $0 [fasta file named with '.fa' and '_' before the fasta to separate the pretss length and the posttss length]\n\n";
die $USAGE unless -e $ARGV[0];
die $USAGE unless $ARGV[0] =~ /.fa$/;
my $fasta = $ARGV[0];
my $name = "$fasta";
$name =~ s/\.fa//;
my @field = split(/\_/,$name);
my $post = splice(@field,-1);
my $pre = splice(@field,-1);
my $out = join('_',@field)."_$pre\_$post\.gc";
open(OUT,">$out");

my $seqio = Bio::SeqIO->new(-file => $fasta,
                            -format => 'fasta');

print OUT "Sequence\tPercentage\tCounts\n";
my @seqs;
while(my $seq = $seqio->next_seq) {
  my($n,$c) = seq_perc($seq->seq);
  print OUT $seq->id."\t".$c."\t".$n."\n";
}

sub seq_perc {

  my $pre_string = shift;
  my $pre_letter = 'CG';

  my $string = uc($pre_string);
  my $letter = uc($pre_letter);
  
  my @n = ($string =~ m/$letter/g);
  my $n = scalar(@n);
  my $n_perc;
 
  if($n) {
    $n_perc = (($n/(length($string)))*100);
  }
  else {
    $n_perc = 0;
  }
  $n_perc =~ s/^(\d+\.\d)\d+$/$1/;
  return($n,$n_perc);
}
