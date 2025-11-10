SAMPLE=$1
CALLER=$2
CALLDIR=$3

../SurVeyor-0.9/bin/compare ../benchmark-data/hgsvc3/$SAMPLE.INSDEL.vcf.gz $CALLDIR/$CALLER/$SAMPLE.vcf.gz -T ../ref/simpleRepeat.hg38.bed -R ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --report --bdup-ids ../benchmark-data/hgsvc3/DUP.ids --force-ids --ignore-seq > $CALLER-ignore_seq/$SAMPLE.txt
