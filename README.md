# surveyor_paper_data

First, unzip and build the bundled SurVeyor. This will prepare the code to benchmark callsets.
```
unzip SurVeyor-0.9.zip 
cd SurVeyor-0.9/
./build_htslib.sh 
cmake -DCMAKE_BUILD_TYPE=Release .
make
```

Download/unzip references
```
cd ref
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
unzip tair10.zip
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.2_ARS-UCD1.3/GCF_002263795.2_ARS-UCD1.3_genomic.fna.gz
gunzip GCF_002263795.2_ARS-UCD1.3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
gunzip -k GCF_000001635.27_GRCm39_genomic.fna.gz
unzip IRGSP-1.0_genomic.fna.zip 
```
