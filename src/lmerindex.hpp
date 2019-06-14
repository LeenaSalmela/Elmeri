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

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <random>

// Parameters
#define MIN_SHARED_MERS 2         // Minimum number of shared mers to consider two Rmaps related
#define QUANTIZE_BAND 1.0         // Width of bin in quantization
#define MIN_FRAGMENT_LEN 1.0      // Minimum fragment length considered
#define MER_SIMILARITY_THRS 5     // Threshold for the similarity of mers to merge their list of Rmaps

#define READ_ID(i) ((i)/2)
#define READ_ORI(i) ((i)%2)

#define FORWARD 0
#define REVERSE 1


// This is needed so we can insert pairs of integers into unordered sets
namespace std {
template <> struct hash<std::pair<int, int>> {
    inline size_t operator()(const std::pair<int, int> &v) const {
        std::hash<int> int_hasher;
        return int_hasher(v.first) ^ int_hasher(v.second);
    }
};
}

// This is needed so we can insert string,integer pairs into unordered sets
namespace std {
  template <> struct hash<std::pair<std::string, int>> {
    inline size_t operator()(const std::pair<std::string, int> &v) const {
      std::hash<std::string> str_hasher;
      std::hash<int> int_hasher;
      return str_hasher(v.first) ^ int_hasher(v.second);
    }
  };
}

// A structure for storing information on an Rmap that is related to a currently inspected Rmap
class related {
public:
  int related_id;    // Id of the related Rmap (i:th Rmap has id 2*i in forward orientation, 2*i+1 in reverse)
  int current_pos;   // Position of the first common mer in the currently inspected Rmap
  int related_pos;   // Position of the first common mer in the related Rmap
};

class lmer_index {
public:
  // Map associating mers (represented as a string) to a set of Rmaps (actually an Rmap index and mer position pair, Rmap i in forward orientation has index 2*i, Rmap i in reverse orientation has index 2*i+1, position is the position of the first fragment in mer)
  std::unordered_map<std::string, std::unordered_set<std::pair<int,int> > *> map;
  // The Rmaps. Rmap i is in position i represented as a vector of fragment lengths
  std::vector<std::vector<double> > rmaps;
  // The names of Rmaps
  std::vector<std::string> names;
  // el or k in case we use k-mers (Note that the minimum number of fragments is the constant MINK when using (el,k)-mers)
  int ell;
  // min k for (el,k)-mer index
  int mink;
  // The gap pattern when using el-mers. Use a gap pattern with no *'s if you want a regular el-mer.
  char *gap_pattern;
  // Threshold for the similarity of mers to merge their list of Rmap
  int mer_similarity_thrs;
  
  // Quantize a fragment length
  int quantize(double d);
  // Extract el-mers (or k-mers) and their starting positions from a sequence of cut sites
  std::unordered_set<std::pair<std::string, int> > extract_lmers(std::vector<double> cuts, int ell, char *gap_pattern);
  // Reads Rmaps from <filename> and creates the index mapping the mers to their Rmaps
  void init(char *filename, int l, int k, char *gp, int sim);
  // Insert a single Rmap into the mer->Rmap map
  void insert(std::vector<double> cuts, int rmap_count);
  // Merge the sets of Rmaps containing a given mer for mers that are similar to each other
  void merge();

  // Get the set of related Rmaps
  int get_related(std::vector<double> rmap, std::vector<std::pair<related, unsigned int> > *counts);
};
