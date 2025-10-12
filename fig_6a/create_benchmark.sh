mkdir truth
cat samples.txt | while read sample ; do
	./create_sample_benchmark.sh $sample
done
