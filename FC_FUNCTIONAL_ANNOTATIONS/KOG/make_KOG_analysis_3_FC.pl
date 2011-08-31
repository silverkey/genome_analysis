#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use Data::Dumper;
use Bio::DB::SeqFeature::Store;

my $user = 'mysql_dev';
my $pwd = 'dEvEl0pEr';
my $GFFdbname = 'FC';
my $WORDdbname = 'FC_upstream_100_1_NSH_100_WL_9';
my $GOfile = 'Fracy1_koginfo_FilteredModels1.tab';
my $word = 'CAACAACAA';
my $outfile = "$word\_KOG_3_$WORDdbname";

my $GFFdb = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                            -user => $user,
                                            -pass => $pwd,
                                            -dsn => "dbi:mysql:$GFFdbname");

my $WORDdb = connect_to_db($WORDdbname,$user,$pwd,'');
my $GO = collect_JGI_GO($GOfile);
my $idmap = get_idmap($GFFdb);
my $resid = get_gene_with_word($WORDdb,$word);
prepare_table($GO,$idmap,$resid,$outfile);

sub prepare_table {
  my $go = shift;
  my $idmap = shift;
  my $resid = shift;
  my $outfile = shift;
  my $total_gene = scalar(keys %{$go->{protein}});
  my $motif_total_gene = get_motif_total_gene_in_go($go,$resid,$idmap);
  open(OUT,">$outfile");
  foreach my $goacc (keys %{$go->{class}}) {
    my $goname = $go->{class}->{$goacc}->{name};
    my $tot_cl = $go->{class}->{$goacc}->{count};
    my $mot_cl = get_motif_gene_in_class($go,$goacc,$resid,$idmap);
    next unless $mot_cl;
    print OUT join("\t","\"$goname\"",$total_gene,$tot_cl,$motif_total_gene,$mot_cl);
    print OUT "\n";
  }
}

sub get_motif_gene_in_class {
  my $go = shift;
  my $goacc = shift;
  my $names = shift;
  my $idmap = shift;
  my $c;
  foreach my $name(@$names) {
    my $id = $idmap->{$name};
    $c ++ if exists $go->{protein}->{$id}->{$goacc};
  }
  return $c;
}

sub get_motif_total_gene_in_go {
  my $go = shift;
  my $names = shift;
  my $idmap = shift;
  my $c;
  foreach my $name(@$names) {
    my $id = $idmap->{$name};
    $c ++ if exists $go->{protein}->{$id};
  }
  return $c;
}

sub collect_JGI_GO {
  my $tab = shift;
  my $href = {};
  open(GO,$tab) or die $!;
  while(my $row = <GO>) {
    chomp($row);
    my($tid,$pid,$kogid,$kogdef,$kogclass,$koggroup) = split(/\t/,$row);
    $href->{protein}->{$pid}->{$koggroup} ++;
    $href->{class}->{$koggroup}->{count} ++;
    $href->{class}->{$koggroup}->{name} = $koggroup;
  }
  return($href);
}

sub get_gene_with_word {
  my $dbh = shift,
  my $word = shift;
  my @id;
  my $id = get_word_id($dbh,$word);
  my $sth = $dbh->prepare('SELECT DISTINCT name FROM seq WHERE id IN (SELECT DISTINCT seq_id FROM seq_occurrence WHERE word_id = '.$id.')');
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    push(@id,$row->{name});
  }
  return \@id;
}

sub get_word_id {
  my $dbh =shift;
  my $word = shift;
  my $sth = $dbh->prepare('SELECT * FROM word WHERE word = '.$dbh->quote($word));
  $sth->execute;
  my $row = $sth->fetchrow_hashref;
  return $row->{id};
}

sub connect_to_db {
  my $db = shift;
  my $usr = shift;
  my $pwd = shift;
  my $host = shift;
  my $dsn = 'dbi:mysql:'.$db;
  $dsn .= ':'.$host if $host; # IN THE CURRENT DBI POD VERSION THERE IS THE '@' IN THE PLACE OF ':'
  my $dbh = DBI->connect($dsn,$usr,$pwd,{PrintWarn=>1,PrintError=>1,RaiseError=>1}) or die $DBI::errstr;
  return $dbh;
}

sub get_idmap {
  my $db = shift;
  my $href = {};
  my @gene = $db->get_features_by_type('gene');
  foreach my $gene(@gene) {
    my %p = $gene->attributes();
    my $alias = $p{Alias}[0];
    $alias =~ s/^t_//;
    my $load_id = $p{load_id}[0];
    $href->{$load_id} = $alias;
  }
  return $href;
}

__END__

Remo-Sangess-MacBook-Pro:GENOME remo$ grep t_135771 Fracy1_GeneModels_FilteredModels1.gff.gff3 
scaffold_5	JGI	gene	1883036	1886348	.	-	.	ID=gene018755;Name=gw1.5.4.1;Alias=t_135771
scaffold_5	JGI	mRNA	1883036	1886348	.	-	.	ID=mRNA018755;Parent=gene018755;Name=t_135771

Remo-Sangess-MacBook-Pro:GENOME remo$ grep 'proteinId 135771' Fracy1_GeneModels_FilteredModels1.gff
scaffold_5	JGI	CDS	1883036	1883886	.	-	0	name "gw1.5.4.1"; proteinId 135771; exonNumber 4
scaffold_5	JGI	CDS	1884008	1884632	.	-	2	name "gw1.5.4.1"; proteinId 135771; exonNumber 3
scaffold_5	JGI	CDS	1884786	1884962	.	-	0	name "gw1.5.4.1"; proteinId 135771; exonNumber 2
scaffold_5	JGI	CDS	1885014	1886348	.	-	0	name "gw1.5.4.1"; proteinId 135771; exonNumber 1

proteinId       gotermId        goName  gotermType      goAcc

# GET gene id from protein id in GFF database from JGI converted GFF
#my @f = $GFFdb->get_features_by_name($alias);
#my @dbid = $f[0]->attributes('parent_id');
#print "$alias -----> $dbid[0]\n";
# GET protein id from gene id in GFF database from JGI converted GFF
#my @f2 = $GFFdb->get_features_by_attribute({load_id => $dbid});
#print Dumper \@f2;
#my @dbid = $f[0]->attributes('parent_id');
#print "$alias -----> $dbid[0]\n";
#print Dumper $GO;
#print "$_\n" foreach $db->types;


t=read.table(file='TATAAA_GO_FC_upstream_100_1_NSH_100_WL_6',sep="\t")
a = c()
for(i in 1:nrow(t)) {
  r = t[i,]
  pval = prop.test(as.numeric(c(r[5],r[3])),as.numeric(c(r[4],r[2])),alternative='g')$p.value
  a = c(a,pval)
}
a = p.adjust(a)
t[a<=0.05,]




t=read.table(file='CAACAA_GO_FC_upstream_100_1_NSH_100_WL_6',sep="\t")
a = c()
for(i in 1:nrow(t)) {
  r = t[i,]
  pval = prop.test(as.numeric(c(r[5],r[3])),as.numeric(c(r[4],r[2])),alternative='g')$p.value
  a = c(a,pval)
}
a = p.adjust(a)
t[a<=0.05,]

