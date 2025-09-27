SAMPLE=$1
CALLER=$2

../SurVeyor-0.9/bin/compare ../benchmark-data/other-orgs/$SAMPLE.vcf.gz ../callsets/ara/$CALLER/$SAMPLE.vcf.gz -T ../ref/TAIR10_chr_all.bed -R ../ref/TAIR10_chr_all.fa --report --bdup-ids ../benchmark-data/other-orgs/$SAMPLE.DUP.ids --force-ids > $CALLER/$SAMPLE.txt
