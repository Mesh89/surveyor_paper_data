CALLER=$1

../SurVeyor-0.9/bin/compare ../benchmark-data/other-orgs/mus_musculus.vcf.gz ../callsets/mus_musculus/$CALLER.vcf.gz -T ../ref/GRCm39_genomic.fna.trf.bed -R ../ref/GCF_000001635.27_GRCm39_genomic.fna --report --bdup-ids ../benchmark-data/other-orgs/mus_musculus.DUP.ids --force-ids > $CALLER/mus_musculus.txt
