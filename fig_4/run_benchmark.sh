for caller in delly manta smoove surveyor surveyor-ara ; do
	mkdir $caller
	cat samples-ara.txt | while read sample ; do
		./compare-ara.sh $sample $caller
	done
done
