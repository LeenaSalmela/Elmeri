/*****************************************************************************
 *  Elmeri
 *
 *  Copyright (C) 2019 Leena Salmela
 *
 *  Contact: leena.salmela@cs.helsinki.fi
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <algorithm>

//#define KMER_INDEX  // Enable this to use a regular k-mer index

#define DOT_OUTPUT // Enable this to output related Rmaps in a dot graph format

#define MINK 5                    // Minimum number of fragments required in an el-mer. Set to 0 for a regular el-mer index.

#include "lmerindex.hpp"

#define MIN(a,b) (((a)<(b))?(a):(b))

/**
 * Command line arguments:
 * - file containing the Rmaps in Valouev format
 * - el (or k)
 * - file for outputting the related Rmaps as a graph in dot format
 * - gap pattern for spaced el-mers (ignored for k-mer index)
 * - threshold for the number of shared mers to consider two Rmaps related
 */
int main(int argc, char *argv[]) {
  if (argc != 7) {
    std::cout << "Invalid arguments." << std::endl;
    exit(2);
  }

  // Command line arguments
  char *rmapfile = argv[1];
  int ell = atoi(argv[2]);
  char *output = argv[3];
  char *gap_pattern = argv[4];
  unsigned int  count_thrs = atoi(argv[5]);
  unsigned int merge_thrs = atoi(argv[6]);
  
  // The el-mer, k-mer or (el,k)-mer index
  lmer_index lind;

  // Read the Rmaps and initialize the index
  lind.init(rmapfile, ell, MINK, gap_pattern, merge_thrs);

#if 0
  // Debugging: print out all Rmaps in the index
  for(int i = 0; i < (int)lind.rmaps.size(); i++) {
    std::cout << i << ":";
    for(auto it = lind.rmaps[i].begin(); it != lind.rmaps[i].end(); ++it) {
      std::cout << " " << (*it);
    }
    std::cout << std::endl;
  }
#endif
  
#if 0
  // Debugging: print out all mers and the indexes of the Rmaps containing them
  for(auto it = lind.map.begin(); it != lind.map.end(); ++it) {
    if (it->second.size() > 1) {
      std::cout << it->first << ":"; 
    
      for(auto it2 = it->second->begin(); it2 != it->second->end(); ++it2) {
	std::cout << " " << it2->first;
      }
      std::cout << std::endl;
    }
  }
#endif

  // Initialize the dot output file
  std::ofstream outfile(output);
#ifdef DOT_OUTPUT
  outfile << "graph {\n";
#endif

  // Gathering statistics on the related Rmaps
  int singletons = 0;
  int edges = 0;
  int nodes = 0;
  
  for(int i = 0; i < (int)lind.rmaps.size(); i++) {
    // Get the related Rmaps
    int counts_size = 0;
    std::vector<std::pair<related, unsigned int> > counts;
    counts_size = lind.get_related(lind.rmaps[i], &counts);
    
    // Output the related Rmaps
    nodes++;
    int degree = 0;

#ifndef DOT_OUTPUT
    outfile << "Rmap " << i << ":";
#endif

    for (int k = 0; k < counts_size; k++) {
      if (counts[k].second >= count_thrs) {
	degree++;
#ifdef DOT_OUTPUT
	// Only print each edge once, i.e. when the current Rmap has a smaller index than the related one
	if ((int)counts[k].first.related_id/2 > i) {
	  // <current Rmap> -- <related Rmap> [weight=<number of shared mers>]
	  outfile << i << " -- " << counts[k].first.related_id/2 << " [weight=" << counts[k].second << "];\n";
	  edges++;
	}
#else
	if ((int)counts[k].first.related_id/2 != i) {
	  // (<related Rmap> <orientation of related Rmap> <pos in related Rmap> <pos in current Rmap>)
	  if (counts[k].first.related_id % 2 == 0) {
	    // Forward Rmap
	    outfile << " (" << counts[k].first.related_id/2 << " F " << counts[k].first.related_pos << " " << counts[k].first.current_pos << ")";
	  } else {
	    // Reverse Rmap
	    outfile << " (" << counts[k].first.related_id/2 << " R " << counts[k].first.related_pos << " " << counts[k].first.current_pos << ")";
	  }
	  edges++;
	}
#endif
      }
    }
#ifndef DOT_OUTPUT
    outfile << std::endl;
#endif

    if (degree == 0) {
      singletons++;
    }
  }

#ifdef DOT_OUTPUT
  outfile << "}\n";
#endif
  outfile.close();

  // Some statistics of the related Rmaps
  printf("Nodes: %d\n", nodes);  // Number of Rmaps
  printf("Edges: %d\n", edges);  // Number of related Rmaps pairs
  printf("Singleton nodes: %d\n", singletons);  // Number of Rmaps with no related Rmaps

  return 0;
}
