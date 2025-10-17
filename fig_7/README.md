Please execute

```
./run_benchmark.sh
python3 ../scripts/generate-xy-benchmark.py . figs/fig figs/callers.txt
```

Fig. 7a will be in figs/

The following instructions will produce the numbers in Fig. 7b

```
mkdir venn
../SurVeyor-0.9/bin/compare nygc.vcf.gz surveyor.vcf.gz --ignore-seq -T ../ref/simpleRepeat.hg38.bed --report -t 32 -f venn/surveyor-only.vcf.gz > venn/nygc-vs-surveyor.txt
../SurVeyor-0.9/bin/compare surveyor.vcf.gz nygc.vcf.gz --ignore-seq -T ../ref/simpleRepeat.hg38.bed --report -t 32 -f venn/nygc-only.vcf.gz > venn/surveyor-vs-nygc.txt
../SurVeyor-0.9/bin/compare variants_GRCh38_sv_insdel_alt_HGSVC2024v1.0.vcf.gz venn/nygc-only.vcf.gz -T ../ref/simpleRepeat.hg38.bed --report -t 32 --ignore-seq > venn/nygc-only.txt
../SurVeyor-0.9/bin/compare variants_GRCh38_sv_insdel_alt_HGSVC2024v1.0.vcf.gz venn/surveyor-only.vcf.gz -T ../ref/simpleRepeat.hg38.bed --report -t 32 --ignore-seq > venn/surveyor-only.txt

echo DEL OVERLAP: `grep "DEL P" venn/nygc-vs-surveyor.txt | tr "/" " " | cut -d" " -f3`
echo matching `grep "DEL P" venn/surveyor-vs-nygc.txt | tr "/" " " | cut -d" " -f3` deletions in NYGC-SV
echo DEL NYGC-SV only: `grep "DEL P" venn/surveyor-vs-nygc.txt | cut -d" " -f3 | awk -F/ '{print $2 - $1}'`
echo validation rate: `cat venn/nygc-only.txt | grep "DEL P" | cut -d" " -f3,4,5`
echo DEL SurVeyor only: `grep "DEL P" venn/nygc-vs-surveyor.txt | cut -d" " -f3 | awk -F/ '{print $2 - $1}'`
echo validation rate: `cat venn/surveyor-only.txt | grep "DEL P" | cut -d" " -f3,4,5`

echo INS+DUP OVERLAP: `grep -E "DUP P|INS P" venn/nygc-vs-surveyor.txt | tr "/" " " | awk '{sum+=$3} END {print sum}'`
echo matching `grep -E "DUP P|INS P" venn/surveyor-vs-nygc.txt | tr "/" " " | awk '{sum+=$3} END {print sum}'` insertions and duplications in NYGC-SV
echo INS+DUP NYGC-SV only: `grep -E "DUP P|INS P" venn/surveyor-vs-nygc.txt | cut -d" " -f3 | awk -F/ '{sum+=$2-$1} END {print sum}'`
echo validation rate: `cat venn/nygc-only.txt | grep -E "DUP P|INS P" | cut -d" " -f3,4,5 | awk -F/ '{n+=$1; d+=$2} END {print n"/"d" = "n/d}'`
echo INS+DUP SurVeyor only: `grep -E "DUP P|INS P" venn/nygc-vs-surveyor.txt | cut -d" " -f3 | awk -F/ '{sum+=$2-$1} END {print sum}'`
echo validation rate: `cat venn/surveyor-only.txt | grep -E "DUP P|INS P" | cut -d" " -f3,4,5 | awk -F/ '{n+=$1; d+=$2} END {print n"/"d" = "n/d}'`
```
