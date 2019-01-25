import sys

bnxfile = sys.argv[1]
enzyme = "BspQI"

j=1
bnx = open(bnxfile, "r");
for line in bnx:
    if line[0] != '#':
        cols = line.split()
        if cols[0] == "1":
            frags = [float(cols[ii]) - float(cols[ii-1]) for ii in range(2,len(cols))]
            if len(frags) >= 10:
                print str(j);
                s = "\t" + enzyme + "\t" + enzyme
                for f in frags:
                    s = s + "\t" + str(f/1000.0)
                print s + "\n"
            j+=1
