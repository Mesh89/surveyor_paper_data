#!/bin/bash

VCF=$1
LR_DATA=$2
WORKDIR=$3
THREADS=$4
OUT_PREFIX=${VCF%.vcf.gz}

# Produce ALT and REF, and use minimap2 to align LR data to them
mkdir -p $WORKDIR
python produce-inversions.py $VCF ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa $WORKDIR/inversions.REF.fa $WORKDIR/inversions.ALT.fa
# if LR data is in BAM or CRAM format, convert it to FASTQ
if [[ $LR_DATA == *.bam || $LR_DATA == *.cram ]]; then
    samtools fastq -@ $THREADS $LR_DATA | gzip > $WORKDIR/long-reads.fq.gz
    LR_DATA=$WORKDIR/long-reads.fq.gz
fi
minimap2 -t $THREADS -ax map-pb -I 48G $WORKDIR/inversions.REF.fa $LR_DATA | samtools view -F 4 -b > $WORKDIR/inversions.REF.bam
minimap2 -t $THREADS -ax map-pb -I 48G $WORKDIR/inversions.ALT.fa $LR_DATA | samtools view -F 4 -b > $WORKDIR/inversions.ALT.bam
rm $WORKDIR/long-reads.fq.gz

# Extract relevant reads, i.e., reads that span a whole REF or ALT contig
python3 extract_spanning_reads.py $WORKDIR/inversions.REF.bam $WORKDIR/inversions.span-full-contig.REF.bam
python3 extract_spanning_reads.py $WORKDIR/inversions.ALT.bam $WORKDIR/inversions.span-full-contig.ALT.bam

python3 dump_fq.py $WORKDIR/inversions.span-full-contig.REF.bam $WORKDIR/inversions.span-full-contig.ALT.bam $WORKDIR/relevant_reads.fq.gz
python3 produce-inversions-naive.py $VCF ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa $WORKDIR/inversions.ALT.large-flanking.fa $WORKDIR/inversions.REF.large-flanking.fa
minimap2 -t $THREADS -ax map-pb -I 48G ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa $WORKDIR/relevant_reads.fq.gz | samtools view -b > $WORKDIR/relevant_reads.FULL-REF.bam
minimap2 -t $THREADS -ax map-pb -I 48G $WORKDIR/inversions.REF.large-flanking.fa $WORKDIR/relevant_reads.fq.gz | samtools view -b > $WORKDIR/relevant_reads.REF.bam
minimap2 -t $THREADS -ax map-pb -I 48G $WORKDIR/inversions.ALT.large-flanking.fa $WORKDIR/relevant_reads.fq.gz | samtools view -b > $WORKDIR/relevant_reads.ALT.bam
python3 classify_vcf.py $WORKDIR/relevant_reads.ALT.bam $WORKDIR/relevant_reads.REF.bam $WORKDIR/relevant_reads.FULL-REF.bam $VCF $OUT_PREFIX.validated-simple-inv.vcf.gz $OUT_PREFIX.validated-with-indel.vcf.gz $OUT_PREFIX.validated-fps.vcf.gz
