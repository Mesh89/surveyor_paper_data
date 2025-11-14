TRUVARI=$1

for caller in delly dragen manta smoove surveyor ; do
	mkdir $caller
	cat samples.txt | while read sample ; do
		./compare-truvari.sh $sample $caller $TRUVARI
	done
done
