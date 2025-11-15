TRUVARI=$1

for caller in surveyor ; do
	cat samples.txt | while read sample ; do
		./compare-truvari.sh $sample $caller ../callsets/hgsvc3-clustered-regt/ $TRUVARI
	done
done

for caller in cutesv sniffles ; do
	cat samples.txt | while read sample ; do
		./compare-truvari.sh $sample $caller ../callsets/hgsvc3/ $TRUVARI
	done
done
