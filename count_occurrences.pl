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
my $out = join('_',@field)."_$pre\_$post\.counts";
open(OUT,">$out");

my $seqio = Bio::SeqIO->new(-file => $fasta,
                            -format => 'fasta');

my @seqs;
while(my $seq = $seqio->next_seq) {
  push(@seqs,$seq);
}

my $counts = nucleotide_count(\@seqs);

print OUT "Position\tA\tT\tC\tG\n";

foreach my $pos (sort {$a <=> $b} keys %$counts) {
  my $rel_pos = $pos-$pre;
  print OUT $rel_pos."\t".$counts->{$pos}->{A}."\t".$counts->{$pos}->{T}."\t".$counts->{$pos}->{C}."\t".$counts->{$pos}->{G}."\n";
}

=head2 nucleotide_count

 Title   : nucleotide_count

 Usage   : my $nuc_count = nucleotide_count($seq_aref);

 Function: count the occurence of nucleotides in a set of sequences

 Returns : hashref with the points in the sequences as keys and the occurrences as values

 Args    : 1 - arrayref of Bio::Seq objects

=cut

sub nucleotide_count {

  my $aref = shift;
  my $counts;

  foreach my $seq(@$aref){
    my @nuc = split('',$seq->seq);
    my $nuc;
    for($nuc=0;$nuc<=$#nuc;$nuc++){
      if($nuc[$nuc] eq 'A' or $nuc eq 'a'){
        $counts->{$nuc}->{'A'}++;
      }
      elsif($nuc[$nuc] eq 'T' or $nuc eq 't'){
        $counts->{$nuc}->{'T'}++;
      }
      elsif($nuc[$nuc] eq 'C' or $nuc eq 'c'){
        $counts->{$nuc}->{'C'}++;
      }
      elsif($nuc[$nuc] eq 'G' or $nuc eq 'g'){
        $counts->{$nuc}->{'G'}++;
      }
    }
  }
  return($counts);
}

