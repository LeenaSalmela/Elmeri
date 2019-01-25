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
