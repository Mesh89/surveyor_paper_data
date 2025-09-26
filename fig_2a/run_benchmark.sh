for caller in delly manta smoove surveyor ; do
	mkdir $caller
	cat samples.txt | while read sample ; do
		./compare.sh $sample $caller
	done
done
