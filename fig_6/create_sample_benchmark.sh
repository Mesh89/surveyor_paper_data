#!/bin/bash

SAMPLE=$1

mkdir truth/$SAMPLE
../SurVeyor-0.9/bin/compare ../benchmark-data/hgsvc3/$SAMPLE.INSDEL.vcf.gz 1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.DEL_DUP.noGT.vcf.gz -T ../ref/simpleRepeat.hg38.bed -R ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --keep-all-called -e | grep -v NONE | grep -v REP | cut -d" " -f1 > truth/$SAMPLE/tps.ids
bcftools view -i "ID==@truth/$SAMPLE/tps.ids" ../benchmark-data/hgsvc3/$SAMPLE.INSDEL.vcf.gz -Oz -o truth/$SAMPLE/tps.vcf.gz
bcftools query -f "%ID\n" -i "SVTYPE=='INS'" truth/$SAMPLE/tps.vcf.gz > truth/$SAMPLE/dup.ids
