for caller in delly manta smoove surveyor ; do
	mkdir $caller
	for depth in 5 10 20 30 40 50 100 ; do
		./compare.sh $depth $caller
	done
done
