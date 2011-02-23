#!/usr/bin/perl

use strict;
use warnings;

use constant CDSEXON => 'CDS';
use constant EXON => 'exon';
use constant MRNA => 'mRNA';
use constant GENE => 'gene';

# MODIFIED FROM GITHUB REPOSITORY OF hyphaltip (AKA Jason Stajich) genome-scripts FOLDER

my $USAGE = "\n\nperl $0 [folder containing jgi-gff files]\n\n";
die $USAGE unless -d $ARGV[0];

my $verbose = 1;
my $gffdir = $ARGV[0];
my $src = 'JGI';

opendir(DIR, $gffdir) || die $!;

FILE: 
for my $file (readdir(DIR)) {
  next unless $file =~ /^(\S+)\.g[tf]f$/;
  my $stem = $1;
  warn("file is $file\n") if $verbose;

  if(-f "$stem\.gff3" && ! -z "$stem\.gff3") {
    warn(" skipping $stem - already processed $stem\.gff3\n");
    next;
  }

  my ($out,$in);
  open($in, "$gffdir/$file") || die $!;
  my (%genes);
  my %lookup;

  while(<$in>) { 

    my ($seqid,$src,$type,$start,$end,$score,$strand,$frame,$lastcol) = split(/\t/,$_,9);
    next unless $type eq 'CDS' || $type eq 'exon';
    
    my ($geneid,$transcriptid);

    if($lastcol =~ /(?:name|gene_id)\s+\"([^\"]+)\";/) {
      $geneid = $1;
    }
    if($lastcol =~ /featureId\s+(\d+)/) {
      $transcriptid = "$1";
    }
    if($lastcol =~ /transcriptId\s+(\d+)/) {
      $transcriptid = "$1";
    } 
    elsif( $lastcol =~ /transcript_id\s+\"([^\"]+)\"/ ) { #"
      $transcriptid = "$1";
    }
	
    if( defined $transcriptid && ! defined $lookup{$geneid} ) {
      $lookup{$geneid} = $transcriptid;
    }
    if( ! defined $transcriptid ) {
      $transcriptid = $lookup{$geneid};
    }
    if( ! defined $geneid || ! defined $transcriptid ) {
      warn("cannot find geneid or transcriptid in $lastcol\n");
      next FILE;
    }

    my $typetag = $type eq 'CDS' ? CDSEXON : EXON;

    push @{$genes{$geneid}->{$transcriptid}}, [$seqid,$typetag,$start,$end,$score,$strand,$frame];
  }

  open( $out, ">$gffdir\/$stem\.gff3" ) || die $!;
  print $out "##gff-version 3\n";
  print $out "#date-prepared ".localtime(time)."\n";

  my %num= ('exon' => 1, 'cds' => 1);

  for my $gene ( keys %genes ) {

    my ($gmin, $gmax,$gstrand,$seqid_all);
    my @lines;

    for my $transcript ( keys %{$genes{$gene}} ) {

      my ($tmin,$tmax,$tstrand);

      for my $CDS ( @{$genes{$gene}->{$transcript}} ) {

        my ($seqid,$type,$start,$end,$score,$strand,$frame) = @$CDS;
        ($start,$end) = sort { $a <=> $b } ($start,$end);
        $gmin = $start if ! defined $gmin || $gmin > $start;
        $gmax = $end if ! defined $gmax || $gmax < $end;
        $gstrand = $strand;

        $tmin = $start if ! defined $tmin || $tmin > $start;
        $tmax = $end if ! defined $tmax || $tmax < $end;
        $tstrand = $strand;

        push @lines, join("\t", $seqid, $src, $type, $start, $end, $score, $strand, $frame,
             join(';',sprintf("ID=%s%05d",lc($type),$num{$type}++),sprintf('Parent=%s',$transcript)));
        $seqid_all = $seqid;
      }

      unshift @lines, join("\t", $seqid_all, $src, MRNA, $tmin, $tmax, '.', $tstrand, '.',
              join(';',sprintf("ID=%s;Name=%s",$transcript,$transcript),sprintf('Parent=%s',$gene)));
    }

    print $out join("\t", $seqid_all, $src, GENE, $gmin, $gmax, '.', $gstrand, '.',
          join(';',sprintf("ID=%s",$gene),sprintf("Name=%s",$gene))),"\n", join("\n", @lines), "\n";
  }
}
