SAMPLE=$1
CALLER=$2
TRUVARI=$3

../scripts/truvari-bench.sh -b ../benchmark-data/hgsvc3/$SAMPLE.INSDEL.vcf.gz -c ../callsets/hgsvc3/$CALLER/$SAMPLE.vcf.gz -w $CALLER/$SAMPLE/ -p $TRUVARI -r ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -d ../benchmark-data/hgsvc3/DUP.ids
cp $CALLER/$SAMPLE/results/summary.txt $CALLER/$SAMPLE.txt
