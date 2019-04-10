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
