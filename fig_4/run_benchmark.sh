for caller in delly manta smoove surveyor surveyor-ara ; do
	mkdir $caller
	cat samples-ara.txt | while read sample ; do
		./compare-ara.sh $sample $caller
	done
done

for caller in delly manta smoove surveyor ; do
	for sample in ERR10310239 ERR10310240 ; do
		./compare-bos_taurus.sh $sample $caller
	done
done

for caller in delly manta smoove surveyor ; do
	./compare-mus_musculus.sh $caller
done

for caller in delly manta smoove surveyor ; do
	for sample in Minghui63 Zhenshan97 ; do
		./compare-oryza_sativa.sh $sample $caller
	done
done

