#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use Bio::Seq;
use DBI;
use Bio::SeqIO;

# This program take in input a fasta file and produce [param] shuffled file of each sequence.
# Than it looks for the word and its tandem repetition in real and shuffled sequences.
# Than it build statistics based on the composition of the real sequences and the shuffled ones

#-----------------------------
# INITIALIZATION
#-----------------------------
my $usage = "\nUsage: perl $0\n
- [fasta file]
- [number of shufflings]
- [word]
- [max number of repetition to check]";
die $usage unless scalar(@ARGV) == 4;
my $fasta = $ARGV[0];
my $nsh = $ARGV[1];
my $word = $ARGV[2];
my $max = $ARGV[3];
# DIRECTORIES WORK
my $startdir = getcwd();
my $workdir = $startdir.'/'.$fasta.'_NSH_'.$nsh.'_'.$word;
$workdir =~ s/\.fa//;
mkdir($workdir) or die "\nmkdir failed: $!\n";
chdir($workdir) or die "\nchdir failed: $!\n";
copy($startdir.'/'.$fasta,$workdir.'/'.$fasta)
or die "\ncopy failed: $!\n";
# OUTPUT TABLES NAME DEFINITION
my $counts_tab = "$fasta\_$nsh\_$word\_counts.xls";
$counts_tab =~ s/\.fasta$//;
$counts_tab =~ s/\.fa$//;
my $fastaname = "$fasta";
$fastaname =~ s/\.fasta$//;
$fastaname =~ s/\.fa$//;

open(OUT,">$counts_tab");
print OUT "rid\tmotif\treal\trandom\n";

my $href = {};

shuffle_fasta($fasta,$nsh);

for(my $c = 1; $c <= $max; $c ++) {
  my $motif = uc($word x $c);
  $href->{$motif}->{real} = '0';
  for(my $x = 1; $x <= $nsh; $x ++) {
    $href->{$motif}->{$x} = '0';
  }
}

my $seqio = Bio::SeqIO->new(-file => $fasta,
                            -format => 'fasta');

while(my $seq = $seqio->next_seq) {
  my $string = uc($seq->seq);
  for(my $c = 1; $c <= $max; $c ++) {
    my $motif = uc($word x $c);
    if($string =~ /$motif/) {
      $href->{$motif}->{real} ++;
    }
  }
}

for(my $x = 1; $x <= $nsh; $x ++) {
  my $sfasta = "$fastaname\_shuffled\_$x\.fa";
  my $seqio = Bio::SeqIO->new(-file => $sfasta,
                              -format => 'fasta');
  while(my $seq = $seqio->next_seq) {
    my $id = $seq->id;
    my $string = uc($seq->seq);
    for(my $c = 1; $c <= $max; $c ++) {
      my $motif = uc($word x $c);
      if($string =~ /$motif/) {
        $href->{$motif}->{$x} ++;
      }
    }
  }
}

my $rn = 1;
foreach my $motif(keys %$href) {
  my $mhref = $href->{$motif};
  my $real = $mhref->{real};
  foreach my $rand(keys %$mhref) {
    next if $rand eq 'real';
    my $random = $mhref->{$rand};
    print OUT join("\t",$rn,$motif,$real,$random);
    print OUT "\n";
    $rn ++;
  }
}

system('rm *.fa');

=head2 shuffle_fasta

Title : shuffle_fasta

Usage : shuffle_fasta($fasta, $shuffling)

Function: create the shuffled versions of fasta

Returns : nothing, the file are created in the working dir
named shuffled_[n].fa where n is the number of the cicle

Args : 1) The fasta to shuffle
2) The number of shifflings

Note : This function use the program "shuffle" from SQUID
the Sean Eddy's C toolkit. Much faster than perl.
http://selab.janelia.org/software.html
The seed used to generate the random numbers are coming
from a counter into the script corresponding to the cicle.
Whatever you need probably Sean Eddy already did....
and it is working better of course ;-)

=cut

sub shuffle_fasta {
  my $fasta = shift;
  my $nsh = shift;
  my $fastaname = "$fasta";
  $fastaname =~ s/\.fasta$//;
  $fastaname =~ s/\.fa$//;
  my $x;
  for($x=1;$x<=$nsh;$x++) {
    my $sfasta = "$fastaname\_shuffled\_$x\.fa";
    system("shuffle --seed $x $fasta > $sfasta");
  }
}
