#!/usr/bin/perl
use strict;
use warnings;
use lib "$ENV{HOME}/src/bioperl-live";
use lib "$ENV{HOME}/src/bioperl-run/lib";
use Bio::Tools::Run::StandAloneBlastPlus;
use Data::Dumper;

$ENV{BLASTPLUSDIR} = "/home/remo/src/ncbi-blast-2.2.25+/bin";

my $file1 = 'Phatr2_chromosomes_geneModels_FilteredModels2_aa.fasta';
my $file2 = 'Thaps3_chromosomes_geneModels_FilteredModels2_aa.fasta';
my $folder = 'PT_vs_TP';
mkdir($folder);
chdir($folder);
system("cp ../$file1 .");
system("cp ../$file2 .");
system("cp ../$0 .");
open(LOG,">LOG") or die "\nERROR opening LOG file: $!\n\n";
system("$ENV{BLASTPLUSDIR}/makeblastdb -in $file1 -dbtype prot -out $file1 -title $file1 -logfile $file1.log");
system("$ENV{BLASTPLUSDIR}/makeblastdb -in $file2 -dbtype prot -out $file2 -title $file2 -logfile $file2.log");

my $out1 = run_blastplus('blastp',$file1,$file2,'blast1');
my $out2 = run_blastplus('blastp',$file2,$file1,'blast2');

my $href1 = parse_blast_results($out1);
my $href2 = parse_blast_results($out2);

open(BRH,">BRH.txt") or die "\nERROR opening BRH file: $!\n\n";

foreach my $id1 (keys %$href1) {
  my $id2 = $href1->{$id1};
  next unless exists $href2->{$id2};
  print BRH "$id1\t$id2\n" if $href2->{$id2} eq $id1;
}

sub run_blastplus {
  my $program = shift;
  my $query = shift;
  my $dbname = shift;
  my $out = shift;

  my $analysis = "$program $query vs $dbname";

  my $exe = $ENV{BLASTPLUSDIR}.'/'.$program;

  my $command = "$exe -query $query -db $dbname -out $out -num_threads 12 -best_hit_overhang 0.1 -evalue 0.01";

  print LOG "\n----\n$command\n";
  print LOG "Executing $analysis\...\n";
  print LOG "Results will be in file $out\n";
  print LOG "$analysis... ";

  system($command);
  sleep(3);

  if(-e $out) { print LOG "OK\n\n" }
  else { print LOG "\nERROR: problems running $analysis\n\n" }
  return $out;
}

sub parse_blast_results {
  my $out = shift;
  my $out_table = "$out\.table";
  my $href = {};
  open(OUT,">$out_table") or die "\nERROR opening $out_table file: $!\n\n";
  print OUT join("\t",qw(clone hit description evalue coverage identity strand frame))."\n";

  my $in = new Bio::SearchIO(-format => 'blast', 
                             -file => $out);

  while(my $result = $in->next_result) {
    my $program = lc($result->algorithm);
    my $candidate = {};
    while(my $hit = $result->next_hit) {
      while(my $hsp = $hit->next_hsp) {
        my $strand = $hsp->strand;
        my $frame = $hsp->frame('query');
        my $identity = $hsp->percent_identity;
        my $coverage = $hsp->length('total') / $result->query_length * 100;
        if($program eq 'blastx') {
          $coverage = $hsp->length('total') / ($result->query_length / 3) * 100;
        }
        $coverage =~ s/^(\d+)\.\d+$/$1/;
        $identity =~ s/^(\d+)\.\d+$/$1/;

        # HARD CODED CUTOFFS!!!
        next unless ($coverage >= 50);

        if((! exists $candidate->{hsp}) || ($hsp->evalue < $candidate->{hsp}->evalue)) {
          $candidate = candidate_this($hit,$hsp,$coverage,$identity,$strand,$frame);
        }
      }
    }

    if(exists $candidate->{hsp}) {
    my $query = strip_id($result->query_name);
    my $accession = strip_id($candidate->{hit}->accession);
    $href->{$query} = $accession;
    print OUT $query ."\t".
              $accession ."\t".
              $candidate->{hit}->description ."\t".
              $candidate->{hsp}->evalue ."\t".
              $candidate->{coverage} ."\t".
              $candidate->{identity}, "\t".
              $candidate->{strand}."\t".
              $candidate->{frame}."\n";
    }
  }
  close(OUT);
  return $href;
}

sub candidate_this {
  my $candidate = {};
  $candidate->{hit} = shift;
  $candidate->{hsp} = shift;
  $candidate->{coverage} = shift;
  $candidate->{identity} = shift;
  $candidate->{strand} = shift;
  $candidate->{frame} = shift;
  return $candidate;
}

sub strip_id {
  my $id = shift;
  my ($acc, $version);
  if ($id =~ /(gb|emb|dbj|sp|tr|pdb|bbs|ref|lcl|tpg)\|(.*)\|(.*)/) {
    ($acc, $version) = split /\./, $2;
  }
  elsif ($id =~ /(pir|prf|pat|gnl)\|(.*)\|(.*)/) {
    ($acc, $version) = split /\./, $3;
  }
  elsif ($id =~ /^jgi\|/) {
    my @f = split(/\|/,$id);
    $acc = $f[2];
  }
  return $acc if $acc;
  return $id;
}
