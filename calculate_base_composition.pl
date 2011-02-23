#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $usage = "\n\tperl $0 [fasta file]\n\n";

die $usage unless -e $ARGV[0];

my $in = $ARGV[0];

my $out = "$in\.composition";
$out =~ s/\.fa//;
open(OUT,">$out");

my $io = Bio::SeqIO->new(-file => $in,
                         -format => 'fasta');

print OUT sprintf("%s\t%s\t%s\t%s\t%s\n",'ID','A','T','C','G');

while(my $seq = $io->next_seq) {
  my $a = seq_perc($seq->seq,'A');
  my $t = seq_perc($seq->seq,'T');
  my $c = seq_perc($seq->seq,'C');
  my $g = seq_perc($seq->seq,'G');
  print OUT sprintf("%s\t%s\t%s\t%s\t%s\n",$seq->id,$a,$t,$c,$g);
}

=head2 seq_perc

 Title   : seq_perc

 Usage   : my $perc = Bio::MCE::Utils->seq_perc($string,$letter);

 Function: it is a function to have the percentage of
           composition of $letter in the $string.
           It take out the N from the calculations

 Returns : a number

 Args    : -1 a string representing a sequence
           -2 a single char representing the letter of wich
              to know the composition

=cut

sub seq_perc {

  my $pre_string = shift;
  my $pre_letter = shift;
  my $N = 'N';

  my $string = uc($pre_string);
  my $letter = uc($pre_letter);
  
  my @n = ($string =~ m/$letter/g);
  my $n = scalar(@n);

  my @nN = ($string =~ m/$N/g);
  my $nN = scalar(@nN);
  
  my $n_perc;
 
  if($n) {
    $n_perc = (($n/(length($string) - $nN)));
  }
  else {
    $n_perc = 0;
  }
  return $n_perc;
}
