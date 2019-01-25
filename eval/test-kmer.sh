# Run elmeri-index for k-mer index (make sure you have compiled elmeri for using these indexes!!!) and evaluate its output


for k in 3 4 5 6 7 8
do
    for th in 2 3 4 5 10
    do
	/usr/bin/time ../src/elmeri-index ../ecoli-sample/ecoli-2000.valouev $k ../ecoli-sample/ecoli-2000-kmer.dot 111 $th
	echo -e "k\tthrs\ttp\tfp\tfn\tprec\trecall"
	echo -n -e "$k\t$th\t"
	python prc.py ../ecoli-sample/ecoli-2000-true.dot ../ecoli-sample/ecoli-2000-kmer.dot
	rm ../ecoli-sample/ecoli-2000-kmer.dot
    done
done

