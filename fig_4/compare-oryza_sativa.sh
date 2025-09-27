SAMPLE=$1
CALLER=$2

../SurVeyor-0.9/bin/compare ../benchmark-data/other-orgs/$SAMPLE.vcf.gz ../callsets/oryza_sativa/$CALLER/$SAMPLE.vcf.gz -T ../ref/IRGSP-1.0_genomic.fna.trf.bed -R ../ref/IRGSP-1.0_genomic.fna --report --bdup-ids ../benchmark-data/other-orgs/$SAMPLE.DUP.ids --force-ids > $CALLER/$SAMPLE.txt
