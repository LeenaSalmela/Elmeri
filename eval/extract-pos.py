# Reads the bed file produced by OMSim and outputs the position in
# fragments and in kbp for each Rmap
# Parameters:
# - Rmaps in Valouev format
# - Bed file
# - reference optical map

import sys

valouevfile=sys.argv[1]
bedfile=sys.argv[2]
reffile=sys.argv[3]

ids = set()

f=open(valouevfile, "r")
i=1
for line in f:
    if i == 1:
        id = int(line)
        ids.add(id)
    i+=1
    if i > 3:
        i=1

f.close()

f=open(reffile, "r")
lines=f.readlines()
cols = lines[1].split()
ref=[]
cumref=[0.0]
for i in range(2,len(cols)):
    ref.append(float(cols[i]))
    cumref.append(float(cols[i]) + cumref[-1])

for i in ref:
    cumref.append(i + cumref[-1])
    
print "rmap start end"
f=open(bedfile, "r")
i=0
for line in f:
    cols = line.split()
    id = int(float(cols[3]))
    if id in ids:
        s = int(cols[1])
        e = int(cols[2])

        sf = -1
        ef = -1
        fsum=0.0
        for j in range(0,len(ref)):
            if sf == -1 and s / 1000.0 < fsum:
                sf = j
            if ef == -1 and e / 1000.0 < fsum:
                ef = j-1
            fsum += ref[j]
        if ef == -1:
            for j in range(0,len(ref)):
                if sf == -1 and s / 1000.0 < fsum:
                    sf = j+len(ref)
                if ef == -1 and e / 1000.0 < fsum:
                    ef = j-1+len(ref)
                fsum += ref[j]

#        print cumref[sf-1], cumref[sf], cumref[ef], cumref[ef+1]
#        print str(i) + " " + str(sf) + " " + str(ef) + " " + str(s/1000.0) + " " + str(e/1000.0)
        print str(i) + " " + str(sf) + " " + str(ef) + " " + str(cumref[sf]) + " " + str(cumref[ef])
        if len(cols[3]) < 3 or cols[3][-2] != ".":
            i+=1
        else:
            if cols[3][-1] == "1":
                i+=1
