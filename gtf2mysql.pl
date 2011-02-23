#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

my $USAGE = "\nperl $0 [gtf file] [dbname]\n\n";

die $USAGE unless -e $ARGV[0];
die $USAGE unless $ARGV[1];

my $gtf_file = $ARGV[0];
my $gff_file = $ARGV[0].'.gff';
my $dbname = $ARGV[1];

my $command1 = 'perl '.$ENV{HOME}.'/src/gtf2gff3/gtf2gff3'." $gtf_file \> $gff_file";
system($command1);

sleep 5;

my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                         -dsn => "dbi:mysql:$dbname",
                                         -user => 'mysql_dev',
                                         -pass => 'riiGbs',
                                         -create => 0);

my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store => $db,
                                                         -verbose => 1,
                                                         -fast => 1,
                                                         -index_subfeatures => 1);

$loader->load($gff_file);
