/*****************************************************************************
 *  Elmeri
 *
 * This file is modified from:
 *
 *  Coral: short reads error correction with multiple alignments
 *  Copyright (C) 2011 Leena Salmela <leena.salmela@cs.helsinki.fi>
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

#include "lmerindex.hpp"
#include "multi-align.hpp"

#include <omp.h>

#include <string.h>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <set>

//#define DEBUG

//#define DOT_OUTPUT

struct sort_pairs {
    bool operator()(const std::pair<int, int> &left, const std::pair<int, int> &right) {
        return left.second > right.second;
    }
};

void usage(char *prog_name) {
    std::cout << "Usage: -i <infile> -o <outfile> -l <ell> -k <min k> -s <space pattern> -c <coverage>" << std::endl;
}

/**
 * Command line arguments:
 * -i <input file> (mandatory)
 * -o <output file> (mandatory)
 * -l <ell> (default: 80)
 * -k <k> (default: 5)
 * -s <spacing pattern> (default: 11111111110001110110010010011101001110001010010100001010011000010111100000001100)
 * -c <coverage> (default: 64, maximum number of Rmaps in an alignment, if set to > 64, will be reset to 64)
 * -t <similarity threshold> (default: 10, threshold for merging similar (el,k)-mers in the index)
 */
int main(int argc, char *argv[]) {

    char *prog_name = argv[0];
    char *infilename = NULL;
    char *outfilename = NULL;

    int ell = 80;
    int mink = 5;
    char *gap_pattern = (char *)"11111111110001110110010010011101001110001010010100001010011000010111100000001100";

    int coverage=ALIGN_SIZE;
    int sim_threshold=MER_SIMILARITY_THRS;
    
    while(argc > 0) {
      if (!strcmp(argv[0], "-i")) {
	if (argc > 1) {
	  infilename = argv[1];
	} else {
	  usage(prog_name);
	  return 1;
	}
	argc--;
	argv++;
      }
      if (!strcmp(argv[0], "-o")) {
	if (argc > 1) {
	  outfilename = argv[1];
	} else {
	  usage(prog_name);
	  return 1;
	}
	argc--;
	argv++;
      }
      if (!strcmp(argv[0], "-l")) {
	if (argc > 1) {
	  ell = atoi(argv[1]);
	} else {
	  usage(prog_name);
	  return 1;
	}
	argc--;
	argv++;
      }
      if (!strcmp(argv[0], "-k")) {
	if (argc > 1) {
	  mink = atoi(argv[1]);
	} else {
	  usage(prog_name);
	  return 1;
	}
	argc--;
	argv++;
      }
      if (!strcmp(argv[0], "-s")) {
	if (argc > 1) {
	  gap_pattern = argv[1];
	} else {
	  usage(prog_name);
	  return 1;
	}
	argc--;
	argv++;
      }
      if (!strcmp(argv[0], "-c")) {
	if (argc > 1) {
	  coverage = atoi(argv[1]);
	  if (coverage > ALIGN_SIZE)
	    coverage = ALIGN_SIZE;
	} else {
	  usage(prog_name);
	  return 1;
	}
	argc--;
	argv++;
      }
      if (!strcmp(argv[0], "-t")) {
	if (argc > 1) {
	  sim_threshold = atoi(argv[1]);
	} else {
	  usage(prog_name);
	  return 1;
	}
	argc--;
	argv++;
      }

      argc--;
      argv++;
    }

    if (infilename == NULL || outfilename == NULL) {
        usage(prog_name);
        exit(1);
    }

    // (el,k)-mer index
    lmer_index lind;

    // Read the Rmaps and initialize the index
    lind.init(infilename, ell, mink, gap_pattern, sim_threshold);
    int total_reads = lind.rmaps.size();

#ifdef DEBUG
    for(auto it = lind.map.begin(); it != lind.map.end(); ++it) {
      // std::string q = it->first;
      // std::cout << q << std::endl;
      std::unordered_set<std::pair<int,int> >* q = it->second;
      std::set<int> q2;
      for(auto it2 = q->begin(); it2 != q->end(); ++it2) {
	q2.insert(it2->first);
      }
      for(auto it3 = q2.begin(); it3 != q2.end(); ++it3) {
	std::cout << "," << *it3;
      }
      std::cout << std::endl;
    }
#endif

    /* Find for each read other reads that shares most (el,k)-mers with it */
    int good_aligns = 0;
    std::ofstream outfile(outfilename);
    int corrected[total_reads];
    for (int i = 0; i < total_reads; i++) {
        corrected[i] = 0;
    }

#ifdef DOT_OUTPUT
    std::cout << "graph {\n";
#endif

    // Set the number of threads
    omp_set_num_threads(NUM_THREADS);

    int current_read = 0;
# pragma omp parallel
    {
      int i = 0;
      while(i < total_reads) {
#pragma omp critical
	{
	  // Find the next read to correct
	  i = current_read;
	  current_read++;
	  if (i % 10000 == 0)
	    std::cout << "Correcting Rmap " << i << std::endl;

	  while(corrected[i] >= 5) {
	    i = current_read;
	    current_read++;
	    if (i % 10000 == 0)
	      std::cout << "Correcting Rmap " << i << std::endl;
	  }
	  
	}
	if (i >= total_reads) {
	  break;
	}
	
	std::vector<double> cuts;
#pragma omp critical
	{
	  cuts.insert(cuts.begin(), lind.rmaps[i].begin(), lind.rmaps[i].end());
	}

	// Get the related Rmaps
	int counts_size = 0;
	std::vector<std::pair<related, unsigned int> > counts;
	related r;
	r.related_id = 2*i;
	r.related_pos = 0;
	r.current_pos = 0;
	// We add the current read here with a very high count to ensure that it will be in the alignment
	counts.push_back(std::make_pair(r, 1000000));

	counts_size = lind.get_related(cuts, &counts);
	
	// Transform the reads to point-patterns (C=cut site in block, A = no cut site in block)
        char *point_patterns[ALIGN_SIZE];
	std::vector<double> cut_vectors[ALIGN_SIZE];
        for (int j = 0; j < counts_size && j < coverage; j++) {
#pragma omp critical
	  {
	    cut_vectors[j].insert(cut_vectors[j].begin(), lind.rmaps[READ_ID(counts[j].first.related_id)].begin(), lind.rmaps[READ_ID(counts[j].first.related_id)].end());
	  }
	}

        int max_len = 0;
        for (int j = 0; j < counts_size && j < coverage; j++) {
#ifdef DEBUG
	  std::cout << "Rmap:";
	  for(int ii = 0; ii < cut_vectors[j].size(); ii++) {
	    std::cout << " " << cut_vectors[j][ii];
	  }
	  std::cout << std::endl;
#endif
	  int rmap_len = cut_vectors[j][cut_vectors[j].size()-1]*1000 + (cut_vectors[j].size()-1)*BLOCK_SIZE;
	  point_patterns[j] = new char[rmap_len / BLOCK_SIZE + 2];
	  if (rmap_len / BLOCK_SIZE + 2 > max_len)
	    max_len = rmap_len / BLOCK_SIZE + 2;
	  int p = 0;
	  int pos = 0;
	  //printf("%d Read: %d Ori: %d\n", friends[j], READ_ID(friends[j]), READ_ORI(friends[j]));
	  if (READ_ORI(counts[j].first.related_id) == FORWARD) {
#ifdef DEBUG
	    std::cout << "F ";
#endif
	    for (int k = 0; k < (int)cut_vectors[j].size(); k++) {
	      if (cut_vectors[j][k]*1000 >= p + BLOCK_SIZE) {
		while (cut_vectors[j][k]*1000 >= p + 2 * BLOCK_SIZE) {
		  point_patterns[j][pos++] = 'A';
		  p += BLOCK_SIZE;
		}
		point_patterns[j][pos++] = 'C';
		p += BLOCK_SIZE;
	      } else {
		point_patterns[j][pos++] = 'C';
	      }
	    }
	    point_patterns[j][pos] = '\0';
	  } else {
#ifdef DEBUG
	    std::cout << "R ";
#endif
	    p = cut_vectors[j][cut_vectors[j].size()-1]*1000;
	    for (int k = cut_vectors[j].size() - 1; k >= 0; k--) {
	      if (cut_vectors[j][k]*1000 <= p - BLOCK_SIZE) {
		while (cut_vectors[j][k]*1000 <= p - 2 * BLOCK_SIZE) {
		  point_patterns[j][pos++] = 'A';
		  p -= BLOCK_SIZE;
		}
		point_patterns[j][pos++] = 'C';
		p -= BLOCK_SIZE;
	      } else {
		point_patterns[j][pos++] = 'C';
	      }
	    }
	    point_patterns[j][pos] = '\0';

	  }
#ifdef DEBUG
	  std::cout << point_patterns[j] << std::endl;
#endif
        }
#ifdef DEBUG
	std::cout << std::endl;
#endif
	
        // Form multiple alignment
        read_pos align_pos[ALIGN_SIZE];
        align *alignment = (align *) malloc(sizeof (align));
        if (alignment == NULL) {
            printf("Out of memory.\n");
            exit(1);
        }
        char *consensus;
        int clen;

        for (int j = 0; j < MIN(coverage, counts_size); j++) {
            align_pos[j].read = j;
            align_pos[j].pos = 0;
            align_pos[j].aligned = 0;
            align_pos[j].edited_read = new char[2 * max_len];
            align_pos[j].edits = 0;
        }

        int num_aligned = multi_align(alignment, point_patterns, align_pos, MIN(coverage, counts_size), 0.5, ALIGN_SIZE, MATCH_REWARD, MM_PENALTY, GAP_PENALTY, BAD_ALIGN_CONSTANT);
        // Realign if some reads (but not too many) where kicked out
        if (num_aligned > 5 && num_aligned < (int)counts.size())
	  num_aligned = multi_align(alignment, point_patterns, align_pos, num_aligned, 0.5, ALIGN_SIZE, MATCH_REWARD, MM_PENALTY, GAP_PENALTY, BAD_ALIGN_CONSTANT);
        clen = get_consensus(alignment, &consensus);
        consensus[clen] = '\0';

#ifdef DEBUG
        // print the alignment
        printf("%i ******************\n", i);

        printf("  %s\n", consensus);

        for (int j = 0; j < num_aligned; j++) {
            if (READ_ORI(counts[align_pos[j].read].first.related_id) == FORWARD)
                printf("F ");
            else
                printf("R ");
            for (int ii = 0; ii < align_pos[j].pos; ii++) {
                printf(" ");
            }
            printf("%s\n", align_pos[j].edited_read);
        }
#endif

        // Check if the alignment is good
        if (num_aligned > GOOD_ALIGN) {
            good_aligns++;
            // good alignment!
            // Count the gaps in the consensus
            int gap_count = 0;
            for (int k = 1; k < clen; k++) {
                if (consensus[k] == 'C')
                    gap_count++;
            }

            // Average gap lengths
            int reads_for_gap[gap_count];
            int avg_gap_len[gap_count];
            int cons_gap_len[gap_count];
            for (int k = 0; k < gap_count; k++) {
                reads_for_gap[k] = 0;
                avg_gap_len[k] = 0;
            }
            int a_count = 0;
            int gap_pos = 0;
            for (int ii = 1; ii < clen; ii++) {
                if (consensus[ii] == 'C') {
                    cons_gap_len[gap_pos++] = (a_count + 1) * BLOCK_SIZE;
                    a_count = 0;
                } else {
                    a_count++;
                }
            }

            for (int k = 0; k < num_aligned; k++) {
                int gap_pos = -1;
                bool prev_cut_ok = false;
                int cut_pos;

                if (READ_ORI(counts[align_pos[k].read].first.related_id) == FORWARD)
                    cut_pos = 0;
                else
  		    cut_pos = cut_vectors[align_pos[k].read].size()-1;

                for (int ii = 0; ii < clen; ii++) {
                    if (consensus[ii] == 'C') {
		      if (align_pos[k].pos <= ii && align_pos[k].pos + (int)strlen(align_pos[k].edited_read) > ii) {
                            if (prev_cut_ok && align_pos[k].edited_read[ii - align_pos[k].pos] == 'C') {
                                reads_for_gap[gap_pos]++;
				if (READ_ORI(counts[align_pos[k].read].first.related_id) == FORWARD)
				  avg_gap_len[gap_pos] += (int) ((cut_vectors[align_pos[k].read][cut_pos] - cut_vectors[align_pos[k].read][cut_pos-1])*1000);
				else
				  avg_gap_len[gap_pos] += (int) ((cut_vectors[align_pos[k].read][cut_pos+1] - cut_vectors[align_pos[k].read][cut_pos])*1000);
                            }
                            if (align_pos[k].edited_read[ii - align_pos[k].pos] == 'C') {
                                prev_cut_ok = true;
                                if (READ_ORI(counts[align_pos[k].read].first.related_id) == FORWARD)
                                    cut_pos++;
                                else
                                    cut_pos--;
                            } else {
                                prev_cut_ok = false;
                            }
                        }
                        gap_pos++;
                    } else {
		      if (align_pos[k].pos <= ii && align_pos[k].pos + (int)strlen(align_pos[k].edited_read) > ii) {
                            if (align_pos[k].edited_read[ii - align_pos[k].pos] == 'C') {
                                prev_cut_ok = false;
                                if (READ_ORI(counts[align_pos[k].read].first.related_id) == FORWARD)
                                    cut_pos++;
                                else
                                    cut_pos--;
                            }
                        }
                    }
                }
            }
	    
            for (int k = 0; k < gap_count; k++) {
                if (reads_for_gap[k] == 0)
                    avg_gap_len[k] = cons_gap_len[k];
                else
                    avg_gap_len[k] = avg_gap_len[k] / reads_for_gap[k];
#ifdef DEBUG
                std::cout << " " << avg_gap_len[k] << "(" << reads_for_gap[k] << ")";
#endif
            }
#ifdef DEBUG
            std::cout << std::endl;
#endif
            // Correct one read.
            for (int j = 0; j < num_aligned; j++) {

#ifdef DOT_OUTPUT
	      // Only print each edge once, i.e. when the current Rmap has a smaller index than the related one
	      int read_id1 = READ_ID(counts[align_pos[0].read].first.related_id);
	      int read_id2 = READ_ID(counts[align_pos[j].read].first.related_id);
	      if (read_id1 < read_id2) {
		// <current Rmap> -- <related Rmap> [weight=<number of shared mers>]
		std::cout << read_id1 << " -- " << read_id2 << " [weight=" << counts[align_pos[j].read].second << "];\n";
	      } else if (read_id2 < read_id1) {
		// <current Rmap> -- <related Rmap> [weight=<number of shared mers>]
		std::cout << read_id2 << " -- " << read_id1 << " [weight=" << counts[align_pos[j].read].second << "];\n";
	      }
	      
#endif

	      
	      //	      if (align_pos[j].read == 0) {
#pragma omp critical
		{
		  int read_id = READ_ID(counts[align_pos[j].read].first.related_id);
#ifdef DEBUG
		  std::cout << read_id << ":";
		  for(auto it = lind.rmaps[read_id].begin(); it != lind.rmaps[read_id].end(); ++it) {
		    std::cout << " " << (*it);
		  }
		  std::cout << std::endl;
#endif
		  lind.rmaps[read_id].clear();
		  lind.rmaps[read_id].push_back(0.0);
		  gap_pos = 0;
		  double cut_pos = 0.0;
		  for (int ii = 1; ii < clen; ii++) {
		    if (consensus[ii] == 'C') {
		      if (align_pos[j].pos < ii && align_pos[j].pos + (int)strlen(align_pos[j].edited_read) > ii) {
			cut_pos += (double) avg_gap_len[gap_pos] / 1000.0;
			lind.rmaps[read_id].push_back(cut_pos);
		      }
		      gap_pos++;
		    }
		  }

		  int ori = READ_ORI(counts[align_pos[j].read].first.related_id);
		  if (ori == REVERSE) {
		    std::vector<double> frags;
		    for (int ii = lind.rmaps[read_id].size()-1; ii > 0; ii--) {
		      frags.push_back(lind.rmaps[read_id][ii] - lind.rmaps[read_id][ii-1]);
		    }
		    lind.rmaps[read_id].clear();
		    lind.rmaps[read_id].push_back(0.0);
		    double cut_pos = 0.0;
		    for(int ii = 0; ii < (int)frags.size(); ii++) {
		      cut_pos += frags[ii];
		      lind.rmaps[read_id].push_back(cut_pos);
		    }
		  }
		  
#ifdef DEBUG
		  std::cout << read_id << ":";
		  for(auto it = lind.rmaps[read_id].begin(); it != lind.rmaps[read_id].end(); ++it) {
		    std::cout << " " << (*it);
		  }
		  std::cout << std::endl;
#endif
		  corrected[READ_ID(counts[align_pos[j].read].first.related_id)] +=1;
		}
		//}
	    }
	}
	    
	for (int j = 0; j < (int)counts.size() && j < coverage; j++) {
	  delete [] align_pos[j].edited_read;
	  delete [] point_patterns[j];
	}
	free(alignment);
      }
    }

    for(int i = 0; i < (int)lind.rmaps.size(); i++) {
        if (corrected[i] > 0)
	  outfile << lind.names[i] << "_C" << std::endl;
        else
	  outfile << lind.names[i] << "_D" << std::endl;
	outfile << "\tEnz\tEnz";
	for (int j = 1; j < (int)lind.rmaps[i].size(); j++) {
	  outfile << "\t" << (lind.rmaps[i][j] - lind.rmaps[i][j-1]);
	}
	outfile << std::endl;
	outfile << std::endl;
    }

#ifdef DOT_OUTPUT
    std::cout << "}\n";
#endif
    
    //std::cout << "Good aligns/Total aligns: " << good_aligns << "/" << reads.size() << std::endl;
#ifndef DOT_OUTPUT
    std::cout << good_aligns << ";" << lind.rmaps.size() << std::endl;
#endif
}
