SAMPLE=$1
CALLER=$2

../SurVeyor-0.9/bin/compare truth/$SAMPLE/tps.vcf.gz ../callsets/nygc-regt/$CALLER/$SAMPLE.vcf.gz -T ../ref/simpleRepeat.hg38.bed -R ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --bdup-ids truth/$SAMPLE/dup.ids --keep-all-benchmark --report --force-ids > $CALLER/$SAMPLE.recall.txt

../SurVeyor-0.9/bin/compare ../benchmark-data/hgsvc3/$SAMPLE.INSDEL.vcf.gz ../callsets/nygc-regt/$CALLER/$SAMPLE.vcf.gz -T ../ref/simpleRepeat.hg38.bed -R ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --report --force-ids > $CALLER/$SAMPLE.precision.txt

grep "DEL SENSITIVITY" $CALLER/$SAMPLE.recall.txt > $CALLER/$SAMPLE.txt
grep "DEL PRECISION" $CALLER/$SAMPLE.precision.txt >> $CALLER/$SAMPLE.txt
grep "DUP SENSITIVITY" $CALLER/$SAMPLE.recall.txt >> $CALLER/$SAMPLE.txt
grep "DUP PRECISION" $CALLER/$SAMPLE.precision.txt >> $CALLER/$SAMPLE.txt

rm $CALLER/$SAMPLE.recall.txt
rm $CALLER/$SAMPLE.precision.txt
