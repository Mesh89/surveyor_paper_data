for caller in delly manta smoove surveyor ; do
	mkdir $caller
	cat samples.txt | while read sample ; do
		./compare.sh $sample $caller
	done
done
exit

for caller in insurveyor+survindel2 ; do
	mkdir $caller
	mkdir surveyor-hgsvc2-shared-samples
	cat hgsvc2-shared-samples.txt | while read sample ; do
		./compare.sh $sample $caller
	done
done
exit

cat hgsvc2-shared-samples.txt | while read sample ; do
	cp surveyor/$sample.txt surveyor-hgsvc2-shared-samples/
        cp surveyor-exclusive/$sample.txt surveyor-hgsvc2-shared-samples-exclusive/
done
