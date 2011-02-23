#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::SeqFeature::Store;

my $dbname = $ARGV[0]; 

my $USAGE = "\nperl $0 [database name]\n\n";

die $USAGE unless $ARGV[0];

my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                         -user => 'mysql_dev',
                                         -pass => 'riiGbs',
                                         -dsn => "dbi:mysql:$dbname");

print "$_\n" foreach $db->types;
