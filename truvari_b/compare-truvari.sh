SAMPLE=$1
CALLER=$2
TRUVARI=$3

../scripts/truvari-bench-compare-seq.sh -b ../benchmark-data/hgsvc3/$SAMPLE.INSDEL.vcf.gz -c ../callsets/hgsvc3-regt/$CALLER/$SAMPLE.vcf.gz -w $CALLER/$SAMPLE/ -r ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -p $TRUVARI -d ../benchmark-data/hgsvc3/DUP.ids -D ../benchmark-data/hgsvc3/DUP.ids
cp $CALLER/$SAMPLE/results/summary.txt $CALLER/$SAMPLE.txt
