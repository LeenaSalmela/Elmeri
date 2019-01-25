# Reads the positions of each Rmap and outputs pairs of Rmaps that
# overlap by at least K fragments and EL kbp. Output is in dot format.
# The Rmaps should be sorted according to their starting position.

import sys

K=7
EL=100.0

i=[]
s=[]
e=[]
sp=[]
ep=[]

# Positions of the Rmaps
# id frag_start frag_end pos_start pos_end
posfile = sys.argv[1]

with open(posfile) as f:
    line = f.readline()
    while line:
        cols = line.split()
        if cols[0] != 'rmap':
            i.append(int(cols[0]))
            s.append(int(cols[1]))
            e.append(int(cols[2]))
            sp.append(float(cols[3]))
            ep.append(float(cols[4]))
        line = f.readline()

posend=4639774
fragend=684

for j in range(0,len(i)):
    for k in range(j+1,len(i)):
        sj = s[j]
        ej = e[j]
        spj = sp[j]
        epj = ep[j]
        sk = s[k]
        ek = e[k]
        spk = sp[k]
        epk = ep[k]

        if epj > posend/1000.0 and epk < spj:
            spk += posend/1000.0
            epk += posend/1000.0
            sk += fragend
            ek += fragend
        if epk > posend/1000.0 and epj < spk:
            spj += posend/1000.0
            epj += posend/1000.0
            sj += fragend
            ej += fragend
 
        if sj < sk:
            if sk < ej:
                if ej < ek:
                    if ej-sk+1 >= K and epj-spk+1 >= EL:
                        if i[j] < i[k]:
                            print str(i[j]) + ' -- ' + str(i[k]) + ';' 
                        else:
                            print str(i[k]) + ' -- ' + str(i[j]) + ';' 
                else:
                    if ek-sk+1 >= K and epk-spk+1 >= EL:
                        if i[j] < i[k]:
                            print str(i[j]) + ' -- ' + str(i[k]) + ';' 
                        else:
                            print str(i[k]) + ' -- ' + str(i[j]) + ';' 
        else:
            if sj < ek:
                if ek < ej:
                    if ek-sj+1 >= K and epk-spj+1 >= EL:
                        if i[j] < i[k]:
                            print str(i[j]) + ' -- ' + str(i[k]) + ';' 
                        else:
                            print str(i[k]) + ' -- ' + str(i[j]) + ';' 
                else:
                    if ej-sj+1 >= K and epj-spj+1 >= EL:
                        if i[j] < i[k]:
                            print str(i[j]) + ' -- ' + str(i[k]) + ';' 
                        else:
                            print str(i[k]) + ' -- ' + str(i[j]) + ';' 
