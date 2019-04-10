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
