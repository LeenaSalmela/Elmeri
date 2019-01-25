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

#include "lmerindex.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <random>

// Quantize a fragment length
int lmer_index::quantize(double d) {
  return (int) (d/QUANTIZE_BAND);
}

#ifndef KMER_INDEX

/**
 *  Extraction of (el,mink)-mers and their starting positions
 */
std::unordered_set<std::pair<std::string, int> > lmer_index::extract_lmers(std::vector<double> cuts, int ell, char *gap_pattern) {
  // The (el,mink)-mers to be returned in the end
  std::unordered_set<std::pair<std::string, int> > lmers;

  // Current offset for the gap pattern in kbp
  double start_pos = 0.0;

  // First cut contained in the current position of the gap pattern
  int first_cut = 0;
  // Last cut contained in the current position of the gap pattern
  int last_cut = 0;

  // First cut that is not masked by the gap pattern
  int first_used_cut = 0;
  
  // Length of the Rmap in kbp
  double len_cuts = cuts[cuts.size()-1];
  // No el-mer can be extracted because the Rmap is too short
  if (len_cuts < ell)
    return lmers;

  // Figure out the last cut contained in the first el-mer
  while (cuts[last_cut] < ell) last_cut++;

  // Cuts left after filtering according to the gap pattern
  std::vector<double> filtered_cuts;
  // Filtered cuts turned into fragment lengths
  std::vector<double> frags;
  // Quantized fragment lengths
  std::vector<int> qfrags;
  // Length (in kbp) of one character in the gap pattern
  double chunk = (double)ell / (double)strlen(gap_pattern);

  // Figure out the first el-mer
  // Filter cuts based on the gap pattern
  first_used_cut = -1;
  for(int jj = first_cut; jj < last_cut; jj++) {
    if (gap_pattern[(int)(cuts[jj]/chunk)] == '1') {
      filtered_cuts.push_back(cuts[jj]);
      if (first_used_cut < 0)
	first_used_cut = jj;
    }
  }
  // Add additional cuts if there are less than mink+1
  for(int jj = last_cut; jj < (int)cuts.size() && filtered_cuts.size() <= (unsigned) mink; jj++) {
    filtered_cuts.push_back(cuts[jj]);
    if (first_used_cut < 0) {
      first_used_cut = jj;
    }
  }
  // Compute the fragment lengths
  for(int jj = 0; jj < (int)filtered_cuts.size()-1; jj++) {
    frags.push_back(filtered_cuts[jj+1]-filtered_cuts[jj]);
  }
  // Quantize the fragment lengths
  for(int jj = 0; jj < (int)frags.size(); jj++) {
    qfrags.push_back(quantize(frags[jj]));
  }
  // Represent the (el,mink)-mer as a string and add it to the list of mers
  std::string qs = "";
  if (qfrags.size() >= (unsigned) mink) {
    for(int jj = 0; jj < (int)qfrags.size(); jj++) {
      qs = qs + "," + std::to_string(qfrags[jj]);
    }
#ifdef DEBUG
    std::cout << "*" << qs << " " << first_used_cut << std::endl;
#endif
    lmers.insert(std::make_pair(qs,first_used_cut));
  }      

  // Add the rest of the (el,mink)-mers
  while (last_cut < (int)cuts.size()) {
    // Find out the minimum shift of the gap pattern that will add or remove a cut site. This can happen at any 1/* or */1 boundary of the gap pattern or at the beginning/end of the gap pattern.

    // Shift needed for the first cut to drop out
    double min_shift = cuts[first_cut]-start_pos;
    // Shift needed for the next cut to enter the gap pattern region
    if (cuts[last_cut] - (start_pos+ell) < min_shift) {
      min_shift = cuts[last_cut] - (start_pos+ell);
    }
    // Minimum shift for any 1/* or */1 boundary to hit the next cut
    for(int jj = 0; jj < (int)strlen(gap_pattern)-1; jj++) {
      if (gap_pattern[jj+1] != gap_pattern[jj]) {
	// A 1/* or */1 boundary
	// -> Find the next cut after that boundary and compute the shift for that cut to hit the boundary
	for(int ii = first_cut; ii < last_cut; ii++) {
	  if (cuts[ii] >= start_pos + (jj+1)*chunk) {
	    if (cuts[ii] - start_pos - (jj+1)*chunk < min_shift) {
	      min_shift = cuts[ii] - start_pos - (jj+1)*chunk;
	      break;
	    }
	  }
	}
      }
    }
#ifdef DEBUG
    std::cout << min_shift << std::endl;
#endif
    // Add a small constant to the shift to make sure that the cut actually changes status
    min_shift += 0.01;
    start_pos += min_shift;

    // Empty the data structures
    frags.clear();
    qfrags.clear();
    filtered_cuts.clear();

    // Find the first cut of this el-mer
    while(first_cut < (int)cuts.size() && start_pos > cuts[first_cut]) first_cut++;
    // Find the last cut of this el-mer
    while(last_cut < (int)cuts.size() && start_pos+ell > cuts[last_cut]) last_cut++;

#ifdef DEBUG
    std::cout << first_cut << " " << last_cut << std::endl;
#endif

    // Check for boundary conditions
    if (first_cut >= (int) cuts.size())
      break;
    if (last_cut >= (int) cuts.size())
      break;

    // Figure out the next el-mer
    // Filter cuts based on the gap pattern
    first_used_cut = -1;
    for(int jj = first_cut; jj < last_cut; jj++) {
      if (gap_pattern[(int)((cuts[jj]-start_pos)/chunk)] == '1') {
	filtered_cuts.push_back(cuts[jj]);
	if (first_used_cut < 0)
	  first_used_cut = jj;
      }
    }
    // Add additional cuts if there are less than mink+1
    for(int jj = last_cut; jj < (int)cuts.size() && filtered_cuts.size() <= (unsigned) mink; jj++) {
      filtered_cuts.push_back(cuts[jj]);
      if (first_used_cut < 0) {
	first_used_cut = jj;
      }
    }
    // Compute the fragment lengths
    for(int jj = 0; jj < (int)filtered_cuts.size()-1; jj++) {
      frags.push_back(filtered_cuts[jj+1]-filtered_cuts[jj]);
    }
    // Quantize the fragment lengths
    for(int jj = 0; jj < (int)frags.size(); jj++) {
      qfrags.push_back(quantize(frags[jj]));
    }
    // Represent the (el,mink)-mer as a string and add it to the list of mers
    std::string qs = "";
    if (qfrags.size() >= (unsigned) mink) {
      for(int jj = 0; jj < (int)qfrags.size(); jj++) {
	qs = qs + "," + std::to_string(qfrags[jj]);
      }
#ifdef DEBUG      
      std::cout << qs << " " << first_used_cut << std::endl;
#endif
      lmers.insert(std::make_pair(qs,first_used_cut));
    }
  }
  
  return lmers;
}
#endif // KMER_INDEX

#ifdef KMER_INDEX

// This version uses regular k-mers (where k==ell and the gap pattern is ignored)
std::unordered_set<std::pair<std::string, int> > lmer_index::extract_lmers(std::vector<double> cuts, int ell, char *gap_pattern) {
  // The k-mers to be returned
  std::unordered_set<std::pair<std::string,int> > lmers;

  // Fragments of the k-mer
  std::vector<double> frags;
  // Quantized fragments of the k-mer
  std::vector<int> qfrags;
    
  for(int i = 0; i < (int)cuts.size()-ell; i++) {
    // Empty the data structures
    frags.clear();
    qfrags.clear();

    // Compute fragment lengths
    for(int jj = 0; jj < ell; jj++) {
      frags.push_back(cuts[i+jj+1]-cuts[i+jj]);
    }
    // Quantize fragment lengths
    for(int jj = 0; jj < (int)frags.size(); jj++) {
      qfrags.push_back(quantize(frags[jj]));
    }
    // Represent the k-mer as a string and add it to the list of k-mers
    std::string qs = "";
    if (qfrags.size() > 1) {
      for(int jj = 0; jj < (int)qfrags.size(); jj++) {
	qs = qs + "," + std::to_string(qfrags[jj]);
      }
#ifdef DEBUG      
      std::cout << qs << std::endl;
#endif
      lmers.insert(std::make_pair(qs,i));
    }
  }

  return lmers;
}
#endif // KMER_INDEX

/**
 * Insert a Rmap represented as a vector of cut sited to the mer map.
 * cuts = the cut sites
 * rma_count = index of the rmap
 */
void lmer_index::insert(std::vector<double> cuts, int rmap_count) {
  // Extract the mers
  std::unordered_set<std::pair<std::string, int> > lmers = extract_lmers(cuts, ell, gap_pattern);

#ifdef DEBUG
  std::cout << lmers.size() << std::endl;
#endif
  
  // Insert the Rmap into the mer map
  for(auto it = lmers.begin(); it != lmers.end(); ++it) {
    std::string qs = it->first;
    int pos = it->second;

#ifdef DEBUG
    std::cout << qs << " " << pos << std::endl;
#endif

    if (map.count(qs) == 1) {
      // A mer that already exists in the map
      map[qs]->insert(std::make_pair(rmap_count,pos));
    } else {
      // A new mer
      std::unordered_set<std::pair<int,int> > *rmaps = new std::unordered_set<std::pair<int,int> >();
      rmaps->insert(std::make_pair(rmap_count, pos));
      map.insert(std::pair<std::string, std::unordered_set<std::pair<int, int> > *>(qs, rmaps));
    }
  }

}

/**
 * Initialize the index
 * - Read Rmaps from filename and insert them to the vector of Rmaps
 * - Construct the mer->Rmap map
 */
void lmer_index::init(char *filename, int l, int k, char *gp) {
  // Save the parameters of the index
  ell = l;
  mink = k;
  gap_pattern = gp;

#ifdef KMER_INDEX
  if (ell < 3 || ell > 10) {
    std::cout << "Invalid value for k: " << ell << std::endl;
    exit(2);
  }
#else
  if (ell < 30 || ell > 200) {
    std::cout << "Invalid value for el: " << ell << std::endl;
    exit(2);
  }
#endif
  
  // Open the input file
  std::ifstream file;
  file.open(filename);
  if (!file.good()) {
    std::cout << "Rmap file invalid." << std::endl;
    exit(1);
  }

  // Read the Rmaps

  std::string line;
  // Line counting: 0=header line (rmap name), 1=fragment line (list of fragments separated by tabs, first fragment is in column 3, 2=empty line)
  int i = 0;
  int rmap_count = 0;
  std::string name;
  
  while(std::getline(file, line)) {
    if (i == 0) {
      name = line;
    }
    if (i == 1) {
      // Read the fragment lengths
      std::vector<std::string> cols;
      boost::split(cols, line, [](char c) {return c == '\t' || c == ' ';});
      if (cols.size() < 3) {
	std::cout << "Invalid fragment list: " << line << std::endl;
	exit(2);
      }
      // Turn the fragment lengths into a sequence of cuts sites
      std::vector<double> cuts;
      double s = 0.0;
      cuts.push_back(s);
      for(int j = 0; j < (int)cols.size(); j++) {
	if (atof(cols[j].c_str()) > 0.0) {
	    s += atof(cols[j].c_str());
	    if (atof(cols[j].c_str()) >= MIN_FRAGMENT_LEN)
	      cuts.push_back(s);
	  }
      }

      // Add the Rmap into the vector of Rmaps
      rmaps.push_back(cuts);
      // Add the name of the Rmap
      names.push_back(name);
      
      // Add the Rmap into the mer->Rmap map
      insert(cuts, rmap_count);
      rmap_count++;

      // Figure out the cut sites for the reverse Rmap
      cuts.clear();
      s = 0.0;
      cuts.push_back(s);
      for(int j = (int)cols.size()-1; j >= 3; j--) {
	s += atof(cols[j].c_str());
	if (atof(cols[j].c_str()) >= MIN_FRAGMENT_LEN)
	  cuts.push_back(s);
      }
      // Add the reversed Rmap into the mer->Rmap map
      insert(cuts, rmap_count);
      rmap_count++;

      if (rmap_count % 10000 == 0)
	std::cout << rmap_count << std::endl;
    }
    
    i++;
    if (i > 2)
      i = 0;
  }

  file.close();

#ifdef DEBUG
  std::cout << map.size() << std::endl;
#endif
  
  // Merge the set of Rmaps for mers that are similar
  merge();
}

#if 0
/**
 * Merge the set of Rmaps for mers that are similar
 * Old version: Makes a copy of the index which is inefficient
 */
void lmer_index::merge() {
  // In each round we merge sets of Rmaps for mers that differ by 1
  int round = 0;
  
  while(round < MER_SIMILARITY_THRS) {
    round++;
    // Use the old sets of Rmaps when merging to avoid propagating longer that distance 1 in each round
    // Therefore we save the old sets here
    std::unordered_map<std::string, std::unordered_set<std::pair<int, int> > > oldmap;
    for(auto it = map.begin(); it != map.end(); ++it) {
      std::string q = it->first;
      std::unordered_set<std::pair<int,int> > old;
      old.insert(it->second->begin(), it->second->end());
      oldmap[q] = old;
    }

    // Iterate over the mers
    for(auto it = oldmap.begin(); it != oldmap.end(); ++it) {
      // The quantized mer
      std::string q = it->first;
      // Split the string to fragment lengths
      std::vector<std::string> cols;
      boost::split(cols, q, [](char c) {return c == ',';});

      // Form all mers such that the quantization of one fragment length is 1 category smaller
      // and merge sets associated to those mers to the set of the current mer
      for(int i=1; i < (int)cols.size(); i++) {
	std::string q2 = "";
	for(int j=1; j < (int)cols.size(); j++) {
	  if (i == j) {
	    int f = atoi(cols[j].c_str());
	    std::stringstream ss;
	    ss << q2 << "," << (f-1);
	    q2 = ss.str();
	  } else {
	    std::stringstream ss;
	    ss << q2 << "," << cols[j];
	    q2 = ss.str();
	  }
	}
	if (oldmap.count(q2) >= 1) {
	  map[q]->insert(oldmap[q2].begin(), oldmap[q2].end());
	}
      }

      // Form all mers such that the quantization of one fragment length is 1 category larger
      // and merge sets associated to those mers to the set of the current mer
      for(int i=1; i < (int)cols.size(); i++) {
	std::string q2 = "";
	for(int j=1; j < (int)cols.size(); j++) {
	  if (i == j) {
	    int f = atoi(cols[j].c_str());
	    std::stringstream ss;
	    ss << q2 << "," << (f+1);
	    q2 = ss.str();
	  } else {
	    std::stringstream ss;
	    ss << q2 << "," << cols[j];
	    q2 = ss.str();
	  }
	}
	if (oldmap.count(q2) >= 1) {
	  map[q]->insert(oldmap[q2].begin(), oldmap[q2].end());
	}
      }
    }
  }
}
#endif

#if 0
/**
 * Merge the set of Rmaps for mers that are similar
 */
void lmer_index::merge() {
  // In each round we merge sets of Rmaps for mers that differ by 1
  int round = 0;

  // A map keeping track which mers have been merged. A -> B means
  // that A has been merged with B. Note that B can have been merged
  // further with some other mer.
  std::unordered_map<std::string, std::string> mergetree;
  
  while(round < MER_SIMILARITY_THRS) {
    round++;
    
    // Iterate over the mers
    for(auto it = map.begin(); it != map.end(); ++it) {
      // The quantized mer
      std::string q = it->first;
      // Split the string to fragment lengths
      std::vector<std::string> cols;
      boost::split(cols, q, [](char c) {return c == ',';});

      // Count the length of the merge path in the mergetree and find the root mer of this mer
      int q_count = 0;
      while(mergetree.count(q) >= 1) {
	q_count++;
	q = mergetree[q];
      }

      // The current mer has already been merged with to a mer that
      // differs from this one by MER_SIMILARITY_THRS. Thus we cannot
      // merge it anymore.
      if (q_count +1 > MER_SIMILARITY_THRS)
	continue;
      
      // Form all mers such that the quantization of one fragment
      // length is 1 category smaller and merge the current mer with
      // those mers if the total merging paths in mergetree will not
      // exceed MER_SIMILARITY_THRS
      for(int i=1; i < (int)cols.size(); i++) {
	std::string q2 = "";
	for(int j=1; j < (int)cols.size(); j++) {
	  if (i == j) {
	    int f = atoi(cols[j].c_str());
	    std::stringstream ss;
	    ss << q2 << "," << (f-1);
	    q2 = ss.str();
	  } else {
	    std::stringstream ss;
	    ss << q2 << "," << cols[j];
	    q2 = ss.str();
	  }
	}
	if (map.count(q2) >= 1) {
	  // Count the length of the merge path in the mergetree and find the root mer of this mer
	  int q2_count = 0;
	  while(mergetree.count(q2) >= 1) {
	    q2_count++;
	    q2 = mergetree[q2];
	  }

	  // Merge the mers if they are not already merged and the
	  // length of the total merge paths does not exceed
	  // MER_SIMILARITY_THRS
	  if (q != q2 && q_count+q2_count+1 <= MER_SIMILARITY_THRS) {
	    // Merge the smaller set to the larger
	    if (map[q2]->size() <= map[q]->size()) {
	      map[q]->insert(map[q2]->begin(), map[q2]->end());
	      mergetree[q2] = q;
	      delete map[q2];
	    } else {
	      map[q2]->insert(map[q]->begin(), map[q]->end());
	      mergetree[q] = q2;
	      delete map[q];
	      q = q2;
	      q_count++;
	      // Check if the current mer now has maximal merge path length 
	      if (q_count+1 > MER_SIMILARITY_THRS)
		break;
	    }
	  }
	}
      }

      // Check if the current mer now has maximal merge path length 
      if (q_count +1 > MER_SIMILARITY_THRS)
	continue;

      // Form all mers such that the quantization of one fragment
      // length is 1 category larger and merge the current mer with
      // those mers if the total merging paths in mergetree will not
      // exceed MER_SIMILARITY_THRS
      for(int i=1; i < (int)cols.size(); i++) {
	std::string q2 = "";
	for(int j=1; j < (int)cols.size(); j++) {
	  if (i == j) {
	    int f = atoi(cols[j].c_str());
	    std::stringstream ss;
	    ss << q2 << "," << (f+1);
	    q2 = ss.str();
	  } else {
	    std::stringstream ss;
	    ss << q2 << "," << cols[j];
	    q2 = ss.str();
	  }
	}
	if (map.count(q2) >= 1) {
	  // Count the length of the merge path in the mergetree and find the root mer of this mer
	  int q2_count = 0;
	  while(mergetree.count(q2) >= 1) {
	    q2_count++;
	    q2 = mergetree[q2];
	  }

	  // Merge the mers if they are not already merged and the
	  // length of the total merge paths does not exceed
	  // MER_SIMILARITY_THRS
	  if (q != q2 && q_count+q2_count+1 <= MER_SIMILARITY_THRS) {
	    // Merge the smaller set to the larger
	    if (map[q2]->size() <= map[q]->size()) {
	      map[q]->insert(map[q2]->begin(), map[q2]->end());
	      mergetree[q2] = q;
	      delete map[q2];
	    } else {
	      map[q2]->insert(map[q]->begin(), map[q]->end());
	      mergetree[q] = q2;
	      delete map[q];
	      q = q2;
	      q_count++;
	      // Check if the current mer now has maximal merge path length 
	      if (q_count+1 > MER_SIMILARITY_THRS)
		break;
	    }
	  }
	}
      }
    }
  }

  // Iterate over the mers and make them point directly to the merged sets
  for(auto it = map.begin(); it != map.end(); ++it) {
      // The quantized mer
      std::string q = it->first;
      if (mergetree.count(q) >= 1) {
	std::string q2 = mergetree[q];
	while(mergetree.count(q2) >= 1)
	  q2 = mergetree[q2];
	map[q] = map[q2];
      }
  }
  
}
#endif

#if 1
/**
 * Merge the set of Rmaps for mers that are similar
 */
void lmer_index::merge() {
  if (MER_SIMILARITY_THRS <= 0)
    return;


  // In each round we merge sets of Rmaps for mers that differ by 1
  int round = 0;
  
  // A map keeping track which mers have been merged. A -> B means
  // that A has been merged with B. Note that B can have been merged
  // further with some other mer.
  std::unordered_map<std::string, std::string> mergetree;
  std::unordered_map<std::string, int> mergesize;

  while(round < 2) {
    round++;
    
    // Iterate over the mers
    for(auto it = map.begin(); it != map.end(); ++it) {
      // The quantized mer
      std::string q = it->first;
      // Split the string to fragment lengths
      std::vector<std::string> cols;
      boost::split(cols, q, [](char c) {return c == ',';});

      // Find the root of the tree
      while(mergetree.count(q) >= 1) {
	q = mergetree[q];
      }
      int qsize = 1;
      if (mergesize.count(q) >= 1) {
	qsize = mergesize[q];
      }

      // The current mer has already been merged too many mers. Thus
      // we cannot merge it anymore.
      if (qsize +1 > MER_SIMILARITY_THRS)
	continue;
      
      // Form all mers such that the quantization of one fragment
      // length is 1 category smaller and merge the current mer with
      // those mers if the total merging paths in mergetree will not
      // exceed MER_SIMILARITY_THRS
      for(int i=1; i < (int)cols.size(); i++) {
	std::string q2 = "";
	for(int j=1; j < (int)cols.size(); j++) {
	  if (i == j) {
	    int f = atoi(cols[j].c_str());
	    std::stringstream ss;
	    ss << q2 << "," << (f-1);
	    q2 = ss.str();
	  } else {
	    std::stringstream ss;
	    ss << q2 << "," << cols[j];
	    q2 = ss.str();
	  }
	}
	if (map.count(q2) >= 1) {
	  // Find the root of the tree
	  while(mergetree.count(q2) >= 1) {
	    q2 = mergetree[q2];
	  }
	  int q2size = 1;
	  if (mergesize.count(q2) >= 1) {
	    q2size = mergesize[q2];
	  }

	  // Merge the mers if they are not already merged and the
	  // size of the merged sets does not exceed MER_SIMILARITY_THRS
	  if (q != q2 && qsize+q2size <= MER_SIMILARITY_THRS) {
	    // Merge the smaller set to the larger
	    if (map[q2]->size() <= map[q]->size()) {
	      map[q]->insert(map[q2]->begin(), map[q2]->end());
	      mergetree[q2] = q;
	      mergesize[q] = qsize+q2size;
	      delete map[q2];
	    } else {
	      map[q2]->insert(map[q]->begin(), map[q]->end());
	      mergetree[q] = q2;
	      mergesize[q2] = qsize+q2size;
	      delete map[q];
	      q = q2;
	    }
	    qsize = qsize+q2size;
	    // Check if the current mer now has maximal merge path length 
	    if (qsize >= MER_SIMILARITY_THRS)
	      break;
	  }
	}
      }

      // Check if the current mer now has maximal merge path length 
      if (qsize >= MER_SIMILARITY_THRS)
	continue;

      // Form all mers such that the quantization of one fragment
      // length is 1 category larger and merge the current mer with
      // those mers if the total merging paths in mergetree will not
      // exceed MER_SIMILARITY_THRS
      for(int i=1; i < (int)cols.size(); i++) {
	std::string q2 = "";
	for(int j=1; j < (int)cols.size(); j++) {
	  if (i == j) {
	    int f = atoi(cols[j].c_str());
	    std::stringstream ss;
	    ss << q2 << "," << (f+1);
	    q2 = ss.str();
	  } else {
	    std::stringstream ss;
	    ss << q2 << "," << cols[j];
	    q2 = ss.str();
	  }
	}
	if (map.count(q2) >= 1) {
	  // Find the root
	  while(mergetree.count(q2) >= 1) {
	    q2 = mergetree[q2];
	  }
	  int q2size = 1;
	  if (mergesize.count(q2) >= 1) {
	    q2size = mergesize[q2];
	  }

	  // Merge the mers if they are not already merged and the
	  // size of the merged sets does not exceed MER_SIMILARITY_THRS
	  if (q != q2 && qsize+q2size <= MER_SIMILARITY_THRS) {
	    // Merge the smaller set to the larger
	    if (map[q2]->size() <= map[q]->size()) {
	      map[q]->insert(map[q2]->begin(), map[q2]->end());
	      mergetree[q2] = q;
	      delete map[q2];
	    } else {
	      map[q2]->insert(map[q]->begin(), map[q]->end());
	      mergetree[q] = q2;
	      delete map[q];
	      q = q2;
	    }
	    qsize = qsize+q2size;
	    // Check if the current mer now has maximal merge path length 
	    if (qsize >= MER_SIMILARITY_THRS)
	      break;
	  }
	}
      }
    }
  }

  // Iterate over the mers and make them point directly to the merged sets
  for(auto it = map.begin(); it != map.end(); ++it) {
      // The quantized mer
      std::string q = it->first;
      if (mergetree.count(q) >= 1) {
	std::string q2 = mergetree[q];
	while(mergetree.count(q2) >= 1)
	  q2 = mergetree[q2];
	map[q] = map[q2];
      }
  }
  
}
#endif


/**
 * Compare two pairs of integers by the second integer of the pairs
 */
struct sort_pairs {
  bool operator()(const std::pair<related, int> &left, const std::pair<related, int> &right) {
    return left.second > right.second;
  }
};

/**
 * Compare two related Rmaps. Sort first by id and secondly by position of first mer in current Rmap
 */
bool sort_related(const related &r1, const related &r2) {
  if (r1.related_id != r2.related_id)
    return r1.related_id < r2.related_id;
  else
    return r1.current_pos < r2.current_pos;
}

/**
 * Get the set of related Rmaps. The related Rmaps are inserted into
 * the counts vector which contains pairs of the relatedness
 * information and the number of el-mers shared by the given Rmap and
 * the related Rmap. Return value is the number of related Rmaps
 * returned.
 */
int lmer_index::get_related(std::vector<double> rmap, std::vector<std::pair<related, unsigned int> > *counts) {
  // Related Rmaps
  std::vector<related> friends;
  // Extract mers for the current Rmap
  std::unordered_set<std::pair<std::string,int> > lmers = extract_lmers(rmap, ell, gap_pattern);

  // If an Rmap has already been added to counts do not readd it
  int exclude_id = counts->size() > 0 ? (*counts)[0].first.related_id : -1;
  
  // Gather the Rmaps related to the current Rmap
  for(auto it = lmers.begin(); it != lmers.end(); ++it) {
    std::string qs = it->first;
    if (map.count(qs) >= 1) {
      // Insert the Rmap ids to friends.
      for(auto it2 = map[qs]->begin(); it2 != map[qs]->end(); it2++) {
	related r;
	r.related_id = it2->first;
	r.current_pos = it->second;
	r.related_pos = it2->second;
	friends.push_back(r);
      }
    }
  }

  // Count how many times each related Rmap is contained in the related Rmaps set
  std::sort(friends.begin(), friends.end(), sort_related);
  related r;
  r.related_id = -1;
  r.related_pos = -1;
  r.current_pos = -1;
  unsigned int c = 0;
  for (int k = 0; k < (int)friends.size(); k++) {
    if (r.related_id >= 0) {
      if (r.related_id != friends[k].related_id) {
	if (r.related_id != exclude_id) {
	  counts->push_back(std::make_pair(r, c));
	}
	r = friends[k];
	c = 1;
      } else {
	c++;
      }
    } else {
      r = friends[k];
      c = 1;
    }
  }
  if (r.related_id != exclude_id) {
    counts->push_back(std::make_pair(r, c));
  }
  
  // Sort the related Rmaps based on the number of shared mers
  std::sort(counts->begin(), counts->end(), sort_pairs());

  return counts->size();
}
