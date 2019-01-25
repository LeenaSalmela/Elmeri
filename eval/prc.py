# Compute precision and recall for output of elmeri
# Usage: python prc.py <dot file of truly related Rmaps> <dot file of predicted related Rmaps>

import re
import sys

# Dot file with truly related Rmaps
tfile=sys.argv[1]
# Dot file of predicted Rmaps
pfile=sys.argv[2]

#tfile='../ecoli/ecoli-sample-true.dot'
#pfile='ecoli-elmeri.dot'

tpairs = dict()
ppairs = dict()

with open(tfile) as f:
    line = f.readline()
    while line:
        m = re.match('(\d+) -- (\d+);', line)
        if m:
            i = int(m.group(1))
            j = int(m.group(2))
            if not i in tpairs:
                tpairs[i] = set()
            tpairs[i].add(j)
        line = f.readline()


with open(pfile) as f:
    line = f.readline()
    while line:
        m = re.match('(\d+) -- (\d+).*;', line)
        if m:
            i = int(m.group(1))
            j = int(m.group(2))
            if not i in ppairs:
                ppairs[i] = set()
            ppairs[i].add(j)
        line = f.readline()

tp = 0
fp = 0
fn = 0

for i in ppairs:
    for j in ppairs[i]:
        if i in tpairs and j in tpairs[i]:
            tp+=1
        else:
            fp+=1

for i in tpairs:
    for j in tpairs[i]:
        if not (i in ppairs and j in ppairs[i]):
            fn += 1


print tp, fp, fn, float(tp)/(float(tp+fp)), float(tp)/(float(tp+fn))
#print 'precision: ' + str(float(tp)/(float(tp+fp)))
#print 'recall: ' + str(float(tp)/(float(tp+fn)))
