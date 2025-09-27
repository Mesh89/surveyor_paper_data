LIBRARY=$1
CALLER=$2

../SurVeyor-0.9/bin/compare ../benchmark-data/other-orgs/bos_taurus.vcf.gz ../callsets/bos_taurus/$CALLER/$LIBRARY.vcf.gz -T ../ref/ARS-UCD1.3_genomic.fna.trf.bed -R ../ref/GCF_002263795.2_ARS-UCD1.3_genomic.fna --report --bdup-ids ../benchmark-data/other-orgs/bos_taurus.DUP.ids --force-ids > $CALLER/$LIBRARY.txt
