#!/usr/bin/perl
use strict;
use warnings;

# GFF3 compliancy of the output validated using
# http://modencode.oicr.on.ca/cgi-bin/validate_gff3_online
#
# Consider to improve the script by using also scoring information

my $input = $ARGV[0];
my $usage = "\nUSGE: perl $0 [newcpgreport output]\n\n";
die $usage unless -e $input;

open(IN,$input);
open(OUT,">$input\.gff3");
print OUT "##gff-version 3\n";

my $id = 0;
my $start = 0;
my $end = 0;

while(my $row = <IN>) {
  my @field = split(/\s+/,$row);
  $id = $field[1] if $field[0] eq 'ID';
  my($start,$end) = split(/\.\./,$field[3]) if $field[0] eq 'FT' && $field[1] eq 'CpG' && $field[2] eq 'island';
  print OUT join("\t",$id,'newcpgreport','CpG_island',$start,$end,'.','.','.','.',"\n") if $start && $end;
  $start = 0;
  $end = 0;
}

__END__

ID   NT_166325  3994 BP.
XX
DE   CpG Island report.
XX
CC   Obs/Exp ratio > 0.60.
CC   % C + % G > 50.00.
CC   Length > 200.
XX
FH   Key              Location/Qualifiers
FT   no islands detected
//
ID   NT_166449  27486 BP.
XX
DE   CpG Island report.
XX
CC   Obs/Exp ratio > 0.60.
CC   % C + % G > 50.00.
CC   Length > 200.
XX
FH   Key              Location/Qualifiers
FT   CpG island       20514..20881
FT                    /size=368
FT                    /Sum C+G=238
FT                    /Percent CG=64.67
FT                    /ObsExp=0.99
FT   CpG island       20978..21328
FT                    /size=351
FT                    /Sum C+G=221
FT                    /Percent CG=62.96
FT                    /ObsExp=0.84
FT   CpG island       21963..22192
FT                    /size=230
FT                    /Sum C+G=152
FT                    /Percent CG=66.09
FT                    /ObsExp=0.79
FT   CpG island       22545..22764
FT                    /size=220
FT                    /Sum C+G=144
FT                    /Percent CG=65.45
FT                    /ObsExp=1.06
FT   CpG island       23221..23705
FT                    /size=485
FT                    /Sum C+G=336
FT                    /Percent CG=69.28
FT                    /ObsExp=0.74
FT   numislands       5
//
