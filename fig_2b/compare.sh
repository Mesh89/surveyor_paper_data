DEPTH=$1
CALLER=$2

../SurVeyor-0.9/bin/compare ../benchmark-data/hgsvc3/NA24385.INSDEL.vcf.gz ../callsets/hg002-giab/$CALLER/$DEPTH.vcf.gz  -T ../ref/simpleRepeat.hg38.bed -R ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --report --bdup-ids ../benchmark-data/hgsvc3/DUP.ids --force-ids > $CALLER/$DEPTH.txt
