TRUVARI=$1

for caller in graphtyper2 paragraph surveyor svtyper ; do
	cat samples.txt | while read sample ; do
		./compare-truvari.sh $sample $caller $TRUVARI
	done
done
