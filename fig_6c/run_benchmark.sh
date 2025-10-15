for caller in graphtyper2 paragraph surveyor svtyper ; do
	mkdir $caller
	cat samples.txt | while read sample ; do
		./compare.sh $sample $caller
	done
done
