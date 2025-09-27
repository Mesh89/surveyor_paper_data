# surveyor_paper_data

First, unzip and build the bundled SurVeyor. This will prepare the code to benchmark callsets.
```
unzip SurVeyor-0.9.zip 
cd SurVeyor-0.9/
./build_htslib.sh 
cmake -DCMAKE_BUILD_TYPE=Release .
make
```

Download the human genome reference and unzip the TAIR10 reference
```
cd ref
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
unzip tair10.zip
```
