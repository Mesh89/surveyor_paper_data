Please execute

```
./run_benchmark.sh
python3 ../scripts/generate-stacked-bar-chart.py figs/samples.txt figs/callers.txt figs/ --separate-legend

mkdir sens-by-SN/
../SurVeyor-0.9/bin/compare ../benchmark-data/hgsvc3/HG00512.INSDEL.vcf.gz ../callsets/hgsvc3/cutesv/HG00512.vcf.gz -T ../ref/simpleRepeat.hg38.bed -R ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --ignore-seq | grep -v NONE | cut -d" " -f1 | sort | uniq > sens-by-SN/cutesv-compare.txt
../SurVeyor-0.9/bin/compare ../benchmark-data/hgsvc3/HG00512.INSDEL.vcf.gz ../callsets/hgsvc3-clustered-regt/surveyor/HG00512.vcf.gz -T ../ref/simpleRepeat.hg38.bed -R ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --ignore-seq | grep -v NONE | cut -d" " -f1 | sort | uniq > sens-by-SN/surveyor-compare.txt
../SurVeyor-0.9/bin/compare ../benchmark-data/hgsvc3/HG00512.INSDEL.vcf.gz ../callsets/hgsvc3/sniffles/HG00512.vcf.gz -T ../ref/simpleRepeat.hg38.bed -R ../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa --ignore-seq | grep -v NONE | cut -d" " -f1 | sort | uniq > sens-by-SN/sniffles-compare.txt
python3 recall-by-sample-number.py ../fig_7/variants_GRCh38_sv_insdel_alt_HGSVC2024v1.0.vcf.gz HG00512 sens-by-SN/cutesv-compare.txt sens-by-SN/sniffles-compare.txt sens-by-SN/surveyor-compare.txt --dup-ids ../benchmark-data/hgsvc3/DUP.ids --outdir sens-by-SN/
```

figs/ will contain Fig. 8a
sens-by-SN/plots/standard fill contain Fig. 8b
