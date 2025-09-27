Please execute

```
mkdir ../callsets/hgsvc3/all-free/
for f in ../callsets/hgsvc3/delly/???????.vcf.gz ; do sample=`basename $f`; echo $sample ; bcftools concat -a ../callsets/hgsvc3/delly/$sample ../callsets/hgsvc3/manta/$sample ../callsets/hgsvc3/smoove/$sample -Oz -o ../callsets/hgsvc3/all-free/$sample ; done 
mkdir ../callsets/hgsvc3/all-free+dragen/
for f in ../callsets/hgsvc3/delly/???????.vcf.gz ; do sample=`basename $f`; echo $sample ; bcftools concat -a ../callsets/hgsvc3/delly/$sample ../callsets/hgsvc3/dragen/$sample ../callsets/hgsvc3/manta/$sample ../callsets/hgsvc3/smoove/$sample -Oz -o ../callsets/hgsvc3/all-free+dragen/$sample ; done
python3 ../scripts/generate-xy-benchmark.py . figs/a/fig figs/a/callers.txt
python3 ../scripts/generate-xy-benchmark.py . figs/b/fig figs/b/callers.txt
```

Figs 3a and 3b will be in figs/
