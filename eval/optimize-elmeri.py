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

# A simulated annealing algorithm to optimize the spacing patterns for Elmeri.
# We assume that:
# - Elmeri is available at ../src/elmeri
# - a set of Rmaps is available at ../ecoli-sample/ecoli-2000.valouev
# - the ground truth of related Rmaps is at ../ecoli-sample/ecoli-2000-true.dot

import os
import random
import math

temp=10000.0
coolingRate=0.003


el=80

seed = ''

for i in range(0,el):
    r = random.randint(0,100)
    if r < 50:
        seed = seed + '0'
    else:
        seed = seed + '1'

os.system('../src/elmeri-index ../ecoli-sample/ecoli-2000.valouev ' + str(el) + " ecoli-2000.dot " + seed + " 2")
os.system('python prc.py ../ecoli-sample/ecoli-2000-true.dot ecoli-2000.dot > prc.txt')

with open('prc.txt') as f:
    lines = f.readlines()
    cols=lines[0].split()
    prc=float(cols[3])
    rec=float(cols[4])

fscore=2*prc*rec/(prc+rec)

print seed, fscore;

optseed=seed
optf=fscore

while(temp > 1):
    r = random.randint(0,len(seed)-1)
    print r
    if seed[r] == '0':
        seed2 = seed[0:r] + '1' + seed[r+1:]
    else:
        seed2 = seed[0:r] + '0' + seed[r+1:]
    os.system('../src/elmeri ../ecoli-sample/ecoli-2000.valouev ' + str(el) + " ecoli-2000.dot " + seed2 + " 2")
    os.system('python prc.py ../ecoli-sample/ecoli-2000-true.dot ecoli-2000.dot > prc.txt')
    with open('prc.txt') as f:
        lines = f.readlines()
        cols=lines[0].split()
        prc2=float(cols[3])
        rec2=float(cols[4])
        fscore2=2*prc2*rec2/(prc2+rec2)
        if fscore2 > optf:
            optseed = seed2
            optf = fscore2
        if fscore2 > fscore:
            seed = seed2
            fscore = fscore2
        else:
            r = random.randint(0,100)
            print fscore2-fscore, r, math.exp(100000*(fscore2-fscore)/temp)*100.0
            if r < math.exp(100000*(fscore2-fscore)/temp)*100.0:
                seed = seed2
                fscore = fscore2
                
    print seed, fscore
    print optseed, optf

    temp = (1.0-coolingRate)*temp
    print temp
