for caller in surveyor ; do
	mkdir $caller
	mkdir $caller-ignore_seq
	cat samples.txt | while read sample ; do
		./compare.sh $sample $caller ../callsets/hgsvc3-clustered-regt/
		./compare-ignore_seq.sh $sample $caller ../callsets/hgsvc3-clustered-regt/
	done
done

for caller in cutesv sniffles ; do
	mkdir $caller
	mkdir $caller-ignore_seq
	cat samples.txt | while read sample ; do
		./compare.sh $sample $caller ../callsets/hgsvc3/
		./compare-ignore_seq.sh $sample $caller ../callsets/hgsvc3/
	done
done
