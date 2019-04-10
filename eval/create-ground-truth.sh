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


# This script reads the positions of Rmaps from a bed file created by OMSim
# and based on them determines the set of related Rmaps

# Reference in silico digested Rmap
REF=../ecoli-sample/MG1655.valouev
# The input Rmaps
RMAPS=../ecoli-sample/ecoli-2000.valouev
# bed file produced by OMSim
BED=../ecoli-sample/ecoli.bed

# Output the positions for each Rmap
python extract-pos.py $RMAPS $BED $REF > ../ecoli-sample/ecoli-2000.pos
# Sort the alignments according to theor starting positions
sort -g -k 2,2 ../ecoli-sample/ecoli-2000.pos > ../ecoli-sample/ecoli-2000.pos.sorted

# Print out the related Rmaps in dot format
python print-ovl.py ../ecoli-sample/ecoli-2000.pos.sorted | sort -k 1,1 -g | awk 'BEGIN{print "graph {";} {print $0;} END {print "}";}' > ../ecoli-sample/ecoli-2000-true.dot
