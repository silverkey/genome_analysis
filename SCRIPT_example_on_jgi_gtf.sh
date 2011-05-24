#!/bin/bash

nohup mysql -umysql_dev -priiGbs < ../scripts/CREATE_DB.mysql > mysql.nohup

nohup perl ../scripts/fa2mysql.pl ../GENOME/Fracy1_assembly_scaffolds_repeatmasked.fasta_N_mask FCm Y > fa2mysql.nohup
nohup perl ../scripts/gtf2mysql.pl ../GENOME/Fracy1_GeneModels_FilteredModels1.gff FCm > gtf2mysql.nohup

nohup perl ../scripts/extract_upstream.pl FCm 100 1 gene > extract_upstream.nohup
nohup perl ../scripts/extract_upstream.pl FCm 200 100 gene > extract_upstream200.nohup
nohup perl ../scripts/extract_upstream.pl FCm 300 200 gene > extract_upstream300.nohup
nohup perl ../scripts/extract_upstream.pl FCm 400 300 gene > extract_upstream400.nohup

nohup perl ../scripts/extract_upstream.pl FCm 30 1 gene > extract_upstream30.nohup
nohup perl ../scripts/extract_upstream.pl FCm 50 20 gene > extract_upstream50.nohup
nohup perl ../scripts/extract_upstream.pl FCm 70 40 gene > extract_upstream70.nohup

nohup perl ../scripts/extract_intergenic.pl FCm gene 1 > extract_intergenic.nohup
nohup perl ../scripts/extract_intronic.pl FCm exon 1 > extract_intronic.nohup
nohup perl ../scripts/extract_exons.pl FCm exon > extract_exon.nohup

nohup perl ../scripts/extract_promoters.pl FCm 1000 500 gene > extract_promoter100.nohup
nohup perl ../scripts/extract_promoters.pl FCm 400 400 gene > extract_promoter400.nohup
nohup perl ../scripts/extract_promoters.pl FCm 200 200 gene > extract_promoter200.nohup
nohup perl ../scripts/extract_promoters.pl FCm 20 20 gene > extract_promoter20.nohup

nohup perl ../scripts/calculate_base_composition.pl FCm_upstream_100_1.fa > composition_upstream.nohup
nohup perl ../scripts/calculate_base_composition.pl FCm_upstream_200_100.fa > composition_upstream2.nohup
nohup perl ../scripts/calculate_base_composition.pl FCm_upstream_300_200.fa > composition_upstream3.nohup
nohup perl ../scripts/calculate_base_composition.pl FCm_upstream_400_300.fa > composition_upstream4.nohup
nohup perl ../scripts/calculate_base_composition.pl FCm_intronic.fa > composition_intronic.nohup
nohup perl ../scripts/calculate_base_composition.pl FCm_intergenic.fa > composition_intergenic.nohup
nohup perl ../scripts/calculate_base_composition.pl FCm_exonic.fa > composition_exonic.nohup

nohup perl ../scripts/count_occurrences.pl FCm_promoter_1000_500.fa > count_occurr_1000.nohup
nohup perl ../scripts/count_gc.pl FCm_promoter_1000_500.fa > count_gc_1000.nohup
nohup perl ../scripts/count_gc.pl FCm_promoter_400_400.fa > count_gc_400.nohup
nohup perl ../scripts/count_gc.pl FCm_promoter_200_200.fa > count_gc_200.nohup
nohup perl ../scripts/count_gc.pl FCm_promoter_20_20.fa > count_gc_20.nohup

R CMD BATCH ../scripts/charts.R
R CMD BATCH ../scripts/charts2.R
