# Run elmeri-index for (el,k)-mer index (make sure you have compiled elmeri for using these indexes!!!) and evaluate its output


#for seed in 11111 11110111111111111011 11100111111111111011 11110111110111111011 11011110110111111011 11011110110111101011 11010110110111101011 11010110110110101101
#for seed in 1101011011011010110111010110110110101101
for seed in 11111111110001110110010010011101001110001010010100001010011000010111100000001100
#for seed in 11010110110110101101
do
    for el in 40 50 60 70 80 90 100
    do
	for th in 2 3 4 5 10
	do
	    /usr/bin/time ../src/elmeri-index ../ecoli-sample/ecoli-2000.valouev $el ../ecoli-sample/ecoli-2000-elmer.dot $seed $th
	    echo -e "el\tthrs\ttp\tfp\tfn\tprec\trecall"
	    echo -n -e "$el\t$th\t"
	    python prc.py ../ecoli-sample/ecoli-2000-true.dot ../ecoli-sample/ecoli-2000-elmer.dot
	    rm ../ecoli-sample/ecoli-2000-elmer.dot
	done
    done
done
