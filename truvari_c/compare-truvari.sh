SAMPLE=$1
CALLER=$2
CALLDIR=$3
TRUVARI=$4

../scripts/truvari-bench.sh -b ../benchmark-data/hgsvc3/$SAMPLE.INSDEL.vcf.gz -c $CALLDIR/$CALLER/$SAMPLE.vcf.gz -r ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa -p $TRUVARI -d ../benchmark-data/hgsvc3/DUP.ids -w $CALLER/$SAMPLE
cp $CALLER/$SAMPLE/results/summary.txt $CALLER/$SAMPLE.txt
