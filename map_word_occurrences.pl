#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;

my $USAGE = "\nUSAGE: perl $0 [fasta file named with '.fa' and '_' before the fasta to separate the pretss length and the posttss length] [word]\n\n";
die $USAGE unless -e $ARGV[0];
die $USAGE unless $ARGV[0] =~ /.fa$/;
die $USAGE unless $ARGV[1] =~ /\w+/;
my $fasta = $ARGV[0];
my $word = uc($ARGV[1]);
my $name = "$fasta";
$name =~ s/\.fa//;
my @field = split(/\_/,$name);
my $post = splice(@field,-1);
my $pre = splice(@field,-1);
my $out = join('_',@field)."_$pre\_$post\_$word\.counts";
open(OUT,">$out");

my $seqio = Bio::SeqIO->new(-file => $fasta,
                            -format => 'fasta');

my @seqs;
while(my $seq = $seqio->next_seq) {
  push(@seqs,$seq);
}

my $counts = word_map(\@seqs,$word);

print OUT "Position\t$word\n";

foreach my $pos (sort {$a <=> $b} keys %$counts) {
  my $rel_pos = $pos-$pre;
  print OUT $rel_pos."\t".$counts->{$pos}."\n";
}

=head2 word_map

 Title   : word_map

 Usage   : my $word_count = ($seq_aref,$word);

 Function: count the occurence of $word in a set of sequences

 Returns : hashref with the points in the sequences as keys and the count as values

 Args    : 1 - arrayref of Bio::Seq objects
           2 - a string representing the word

=cut

sub word_map {
  my $aref = shift;
  my $word = shift;
  my $counts;

  foreach my $seq(@$aref){
    my $string = uc($seq->seq);
    my $offset = 0;
    my $result = index($string,$word,$offset);

    while ($result != -1) {
      $counts->{$result} ++;
      $offset = $result + 1;
      $result = index($string,$word,$offset);
    }
  }
  return($counts);
}

__END__
t=read.table(file='FC_M1_promoter_1000_500_CAA.counts',header=T)
plot(t[,1],t[,2],type='l')
