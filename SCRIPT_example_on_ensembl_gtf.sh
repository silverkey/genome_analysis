#!/bin/bash

nohup mysql -umysql_dev -priiGbs < ../scripts/CREATE_DB.mysql > mysql.nohup

nohup perl ../scripts/fa2mysql.pl ../GENOME/Thalassiosira_pseudonana.Thaps3.dna_rm.toplevel.fa TPm Y > fa2mysql.nohup
nohup perl ../scripts/gtf2mysql.pl ../GENOME/thalassiosira_pseudonana.Thaps3.59.gtf TPm > gtf2mysql.nohup

nohup perl ../scripts/extract_upstream.pl TPm 100 1 gene:protein_coding > extract_upstream.nohup
nohup perl ../scripts/extract_upstream.pl TPm 200 100 gene:protein_coding > extract_upstream200.nohup
nohup perl ../scripts/extract_upstream.pl TPm 300 200 gene:protein_coding > extract_upstream300.nohup
nohup perl ../scripts/extract_upstream.pl TPm 400 300 gene:protein_coding > extract_upstream400.nohup

nohup perl ../scripts/extract_upstream.pl TPm 30 1 gene:protein_coding > extract_upstream30.nohup
nohup perl ../scripts/extract_upstream.pl TPm 50 20 gene:protein_coding > extract_upstream50.nohup
nohup perl ../scripts/extract_upstream.pl TPm 70 40 gene:protein_coding > extract_upstream70.nohup

nohup perl ../scripts/extract_intergenic.pl TPm gene:protein_coding 1 > extract_intergenic.nohup
nohup perl ../scripts/extract_intronic.pl TPm exon:protein_coding 1 > extract_intronic.nohup
nohup perl ../scripts/extract_exons.pl TPm exon:protein_coding > extract_exon.nohup

nohup perl ../scripts/extract_promoters.pl TPm 1000 500 gene:protein_coding > extract_promoter100.nohup
nohup perl ../scripts/extract_promoters.pl TPm 400 400 gene:protein_coding > extract_promoter400.nohup
nohup perl ../scripts/extract_promoters.pl TPm 200 200 gene:protein_coding > extract_promoter200.nohup
nohup perl ../scripts/extract_promoters.pl TPm 20 20 gene:protein_coding > extract_promoter20.nohup

nohup perl ../scripts/calculate_base_composition.pl TPm_upstream_100_1.fa > composition_upstream.nohup
nohup perl ../scripts/calculate_base_composition.pl TPm_upstream_200_100.fa > composition_upstream2.nohup
nohup perl ../scripts/calculate_base_composition.pl TPm_upstream_300_200.fa > composition_upstream3.nohup
nohup perl ../scripts/calculate_base_composition.pl TPm_upstream_400_300.fa > composition_upstream4.nohup
nohup perl ../scripts/calculate_base_composition.pl TPm_intronic.fa > composition_intronic.nohup
nohup perl ../scripts/calculate_base_composition.pl TPm_intergenic.fa > composition_intergenic.nohup
nohup perl ../scripts/calculate_base_composition.pl TPm_exonic.fa > composition_exonic.nohup

nohup perl ../scripts/count_occurrences.pl TPm_promoter_1000_500.fa > count_occurr_1000.nohup
nohup perl ../scripts/count_gc.pl TPm_promoter_1000_500.fa > count_gc_1000.nohup
nohup perl ../scripts/count_gc.pl TPm_promoter_400_400.fa > count_gc_400.nohup
nohup perl ../scripts/count_gc.pl TPm_promoter_200_200.fa > count_gc_200.nohup
nohup perl ../scripts/count_gc.pl TPm_promoter_20_20.fa > count_gc_20.nohup

R CMD BATCH ../scripts/charts.R
R CMD BATCH ../scripts/charts2.R
