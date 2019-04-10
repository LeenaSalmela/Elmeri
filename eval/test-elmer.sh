#  Elmeri
#
#  Copyright (C) 2019 Leena Salmela
#
#  Contact: leena.salmela@cs.helsinki.fi
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


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
