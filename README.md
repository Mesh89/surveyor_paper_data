# surveyor_paper_data

First, unzip and build the bundled SurVeyor. This will prepare the code to benchmark callsets.
```
unzip SurVeyor-0.9.zip 
cd SurVeyor-0.9/
./build_htslib.sh 
cmake -DCMAKE_BUILD_TYPE=Release .
make
```