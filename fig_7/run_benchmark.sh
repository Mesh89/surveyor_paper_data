for caller in nygc surveyor ; do
	mkdir $caller
	cat samples.txt | while read sample ; do
		./compare-ignore_seq.sh $sample $caller
	done
done
