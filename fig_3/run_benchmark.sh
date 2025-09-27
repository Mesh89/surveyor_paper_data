for caller in dragen surveyor all-free all-free+dragen ; do
	mkdir $caller
	cat samples.txt | while read sample ; do
		./compare.sh $sample $caller
	done
done
