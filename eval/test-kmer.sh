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

