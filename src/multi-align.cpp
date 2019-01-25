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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "multi-align.hpp"

//#define DEBUG_ALIGN

/* Update the consensus sequence of an alignment starting at position
   start and ending at position end */
void update_consensus(align *alignment, int start, int end) {
    int i;

    if (end >= MAX_CONTIG_LEN) {
        printf("multi-align: max contig len too short.\n");
        exit(1);
    }

#ifdef DEBUG_ALIGN
    printf("Updating consensus: %d -> %d\n", start, end);
#endif

    for (i = start; i < end; i++) {
        if ((alignment->contigA[i] == 0 && alignment->contigC[i] == 0 &&
                alignment->contigG[i] == 0 && alignment->contigT[i] == 0) ||
                (alignment->contigN[i] > alignment->contigA[i] &&
                alignment->contigN[i] > alignment->contigC[i] &&
                alignment->contigN[i] > alignment->contigG[i] &&
                alignment->contigN[i] > alignment->contigT[i])) {
            alignment->consensus[i] = '-';
        } else if (alignment->contigA[i] > alignment->contigC[i] &&
                alignment->contigA[i] > alignment->contigG[i] &&
                alignment->contigA[i] > alignment->contigT[i]) {
            alignment->consensus[i] = 'A';
        } else if (alignment->contigC[i] > alignment->contigG[i] &&
                alignment->contigC[i] > alignment->contigT[i]) {
            alignment->consensus[i] = 'C';
        } else if (alignment->contigG[i] > alignment->contigT[i]) {
            alignment->consensus[i] = 'G';
        } else {
            alignment->consensus[i] = 'T';
        }
    }
}

/**
 * Align two reads with the given scoring scheme and return the number
 * of added/deleted C's in the alignment. Scoring scheme includes the
 * weight alph_skew which is a multiplier for match_reward for C's
 */
int align_two_reads(align *alignment, char **reads, int i1, int i2,
		    int match_reward, int mm_penalty, int gap_penalty, int alph_skew) {
  // Pointers to the reads
  char *r1 = reads[i1];
  char *r2 = reads[i2];

  // Read lengths
  int m1 = strlen(r1);
  int m2 = strlen(r2);

  int i,j;

#ifdef DEBUG_ALIGN
  printf("%s\n", r1);
  printf("%s\n", r2);
#endif

  // Initialize dynamic programming array: free edits in the beginning
  // and end of the alignments (i.e. do an overlap alignment)
  for (i = 0; i <= m1; i++) {
    alignment->dp[i][0] = 0;
    alignment->dp_trace[i][0] = 'J';
  }
  for (j = 0; j <= m2; j++) {
    alignment->dp[0][j] = 0;
    alignment->dp_trace[0][j] = 'I';
  }

  // Fill in the dp matrix
  for (i = 0; i < m1; i++) {
    // Gaps are free if one of the sequences has no more chars left
    int gapi = (i == m1 - 1) ? 0 : gap_penalty;
    for (j = 0; j < m2; j++) {
      // Gaps are free if one of the sequences has no more chars left
      int gapj = (j == m2 - 1) ? 0 : gap_penalty;
      if (toupper(r2[j]) == 'A') {
	int mm = (toupper(r1[i]) == toupper(r2[j])) ? match_reward : -mm_penalty;
	if (alignment->dp[i][j] + mm >= alignment->dp[i][j + 1] - gapj &&
	    alignment->dp[i][j] + mm >= alignment->dp[i + 1][j] - gapi) {
	  alignment->dp[i + 1][j + 1] = alignment->dp[i][j] + mm;
	  if (r1[i] == 'A') {
	    alignment->dp_trace[i + 1][j + 1] = 'M';
	  } else {
	    alignment->dp_trace[i + 1][j + 1] = 'D';
	  }
	} else if (alignment->dp[i][j + 1] - gapj >= alignment->dp[i + 1][j] - gapi) {
	  alignment->dp[i + 1][j + 1] = alignment->dp[i][j + 1] - gapj;
	  alignment->dp_trace[i + 1][j + 1] = 'J';
	} else {
	  alignment->dp[i + 1][j + 1] = alignment->dp[i + 1][j] - gapi;
	  alignment->dp_trace[i + 1][j + 1] = 'I';
	}
      } else {
	int mm = (toupper(r1[i]) == toupper(r2[j])) ? alph_skew*match_reward : -mm_penalty;
	if (alignment->dp[i][j] + mm >= alignment->dp[i][j + 1] - gapj &&
	    alignment->dp[i][j] + mm >= alignment->dp[i + 1][j] - gapi) {
	  alignment->dp[i + 1][j + 1] = alignment->dp[i][j] + mm;
	  if (r1[i] == 'C') {
	    alignment->dp_trace[i + 1][j + 1] = 'M';
	  } else {
	    alignment->dp_trace[i + 1][j + 1] = 'D';
	  }
	} else if (alignment->dp[i][j + 1] - gapj >= alignment->dp[i + 1][j] - gapi) {
	  alignment->dp[i + 1][j + 1] = alignment->dp[i][j + 1] - gapj;
	  alignment->dp_trace[i + 1][j + 1] = 'J';
	} else {
	  alignment->dp[i + 1][j + 1] = alignment->dp[i + 1][j] - gapi;
	  alignment->dp_trace[i + 1][j + 1] = 'I';
	}
      }
    }
  }

  /* Traceback in the dp array to count edits */
  i = m1;
  j = m2;
  int edits = 0;

  while (j > 0) {
    if (alignment->dp_trace[i][j] == 'D' ||
	alignment->dp_trace[i][j] == 'M') {
#ifdef DEBUG_ALIGN
      printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
      // A mismatch always involves a C
      if (alignment->dp_trace[i][j] == 'D') {
	edits++;
      }
      i--;
      j--;
    } else if (alignment->dp_trace[i][j] == 'J') {
#ifdef DEBUG_ALIGN
      printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
      // Overlap alignment -> edits in the beginning/end are free
      if (j > 0 && j < m2) {
	if (r1[i-1] == 'C')
	  edits++;
      }
      i--;
    } else if (alignment->dp_trace[i][j] == 'I') {
#ifdef DEBUG_ALIGN
      printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
      // Overlap alignment -> edits in the beginning/end are free
      if (i > 0 && i < m1) {
	if (r2[j-1] == 'C')
	  edits++;
      }
      j--;
    } else {
      printf("Undefined dp_trace %d %d\n", i, j);
    }
  }

#ifdef DEBUG_ALIGN
  printf("\nTracing done. %d\n", edits);
#endif
  
  return edits;
}



/**
 * Align one more read against the alignment.
 * - alignment: The alignment computed so far
 * - reads: The set of reads
 * - pos: The set of reads for which the alignment is computed
 * - end: Length of alignment
 * - p: The id to pos of the new read to align (the reads with smaller ids have already been aligned)
 * - max_error_rate: if more errors, set the read as unaligned
 * - scoring scheme: match of A: match_reward, match of C: alph_skew*match_reward, 
 *                   mismatch cost: mm_penalty, indel cost: gap_penalty
 * Returns the length of the alignment after the new read has been added.
 */
int align_read(align *alignment, char **reads, read_pos *pos,
        int end, int p, double max_error_rate,
        int match_reward, int mm_penalty, int gap_penalty, int alph_skew) {
    int i, j;

    // Read and read length
    char *read = reads[pos[p].read];
    int m = strlen(read);

    // Buffer for edited read and its length
    char er[MAX_READ_LENGTH];
    int len;

    int k;
    int kk;
    int m2;

    int edits;
    int old_end;

    if (end <= 0) {
        // length of consensus is 0
        // update base counts
        if (m >= MAX_CONTIG_LEN) {
            printf("Max contig length too small\n");
            exit(1);
        }
        for (i = 0; i < (signed int) m; i++) {

            alignment->contigA[i] = 0;
            alignment->contigC[i] = 0;
            alignment->contigG[i] = 0;
            alignment->contigT[i] = 0;
            alignment->contigN[i] = 0;
            switch (read[i]) {
                case 'a':
                case 'A':
                    alignment->contigA[i]++;
                    break;
                case 'c':
                case 'C':
                    alignment->contigC[i]++;
                    break;
                case 'g':
                case 'G':
                    alignment->contigG[i]++;
                    break;
                case 't':
                case 'T':
                    alignment->contigT[i]++;
                    break;
                case '-':
                case'N':
                case 'n':
                    alignment->contigN[i]++;
                    break;
            }
        }
        alignment->offset = pos[p].pos;
        pos[p].pos = 0;
        strcpy(pos[p].edited_read, reads[pos[p].read]);
        pos[p].aligned = 1;
        return m;
    }

    update_consensus(alignment, 0, end);
    alignment->consensus[end] = '\0';

    pos[p].pos = pos[p].pos - alignment->offset;


#ifdef DEBUG_ALIGN
    printf("Aligning read %s (%d, %d)\n", read, pos[p].read, pos[p].pos);
    printf("Against conse %s\n", alignment->consensus);
#endif

    // Full alignment

    /* Align the read against the consensus */

    // Initialize dynamic programming matrix
    // Overlap alignment -> Free indels in the beginning of the sequences
    for (i = 0; i <= end; i++) {
        alignment->dp[i][0] = 0;
        alignment->dp_trace[i][0] = 'J';
    }
    for (j = 0; j <= m; j++) {
        alignment->dp[0][j] = 0;
        alignment->dp_trace[0][j] = 'I';
    }

    // Fill the dp matrix
    for (i = 0; i < end; i++) {
        // Insertions are free if the consensus reaches no further
        int gapi = (i == end - 1) ? 0 : (toupper(alignment->consensus[i]) == 'A' ? gap_penalty : gap_penalty);
        for (j = 0; j < m; j++) {
            // Insertions are free if the read reaches no further
            int gapj = (j == m - 1 || alignment->consensus[i] == '-') ? 0 : (toupper(read[j]) == 'A' ? gap_penalty : gap_penalty);
            if (toupper(read[j]) == 'A') {
                int mm = (alignment->contigA[i] * match_reward - alignment->contigC[i] * mm_penalty - alignment->contigN[i] * mm_penalty) / (alignment->contigA[i] + alignment->contigC[i] + alignment->contigN[i]);
                if (alignment->dp[i][j] + mm >= alignment->dp[i][j + 1] - gapj &&
                        alignment->dp[i][j] + mm >= alignment->dp[i + 1][j] - gapi) {
                    alignment->dp[i + 1][j + 1] = alignment->dp[i][j] + mm;
                    if (alignment->consensus[i] == 'A') {
                        alignment->dp_trace[i + 1][j + 1] = 'M';
                    } else {
                        alignment->dp_trace[i + 1][j + 1] = 'D';
                    }
                } else if (alignment->dp[i][j + 1] - gapj >= alignment->dp[i + 1][j] - gapi) {
                    alignment->dp[i + 1][j + 1] = alignment->dp[i][j + 1] - gapj;
                    alignment->dp_trace[i + 1][j + 1] = 'J';
                } else {
                    alignment->dp[i + 1][j + 1] = alignment->dp[i + 1][j] - gapi;
                    alignment->dp_trace[i + 1][j + 1] = 'I';
                }
            } else {
                int mm = (alph_skew * alignment->contigC[i] * match_reward - alignment->contigA[i] * mm_penalty - alignment->contigN[i] * mm_penalty) / (alignment->contigA[i] + alignment->contigC[i] + alignment->contigN[i]);
                if (alignment->dp[i][j] + mm >= alignment->dp[i][j + 1] - gapj &&
                        alignment->dp[i][j] + mm >= alignment->dp[i + 1][j] - gapi) {
                    alignment->dp[i + 1][j + 1] = alignment->dp[i][j] + mm;
                    if (alignment->consensus[i] == 'C') {
                        alignment->dp_trace[i + 1][j + 1] = 'M';
                    } else {
                        alignment->dp_trace[i + 1][j + 1] = 'D';
                    }
                } else if (alignment->dp[i][j + 1] - gapj >= alignment->dp[i + 1][j] - gapi) {
                    alignment->dp[i + 1][j + 1] = alignment->dp[i][j + 1] - gapj;
                    alignment->dp_trace[i + 1][j + 1] = 'J';
                } else {
                    alignment->dp[i + 1][j + 1] = alignment->dp[i + 1][j] - gapi;
                    alignment->dp_trace[i + 1][j + 1] = 'I';
                }
            }
        }
    }

    /* Traceback in the dp array - first count edits */
    i = end;
    j = m;
    len = 0;
    edits = 0;

    while (i >= 0 && alignment->dp_trace[i][j] == 'J') {
#ifdef DEBUG_ALIGN
        printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
        i--;
    }

    int insC=0;
    int delC=0;
    while (j > 0) {
        if (alignment->dp_trace[i][j] == 'D' ||
                alignment->dp_trace[i][j] == 'M') {
#ifdef DEBUG_ALIGN
            printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
            len++;
            if (alignment->dp_trace[i][j] == 'D') {
                edits++;
		if (read[j-1] == 'C')
		  insC++;
		else if (alignment->consensus[i-1] == 'C')
		  delC++;
	    }
            i--;
            j--;
        } else if (alignment->dp_trace[i][j] == 'J') {
#ifdef DEBUG_ALIGN
            printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
            if (j > 0 && j < m && alignment->consensus[i - 1] != '-') {
                len++;
                edits++;
		if (alignment->consensus[i-1] == 'C')
		  delC++;
            }
            i--;
        } else if (alignment->dp_trace[i][j] == 'I') {
#ifdef DEBUG_ALIGN
            printf("%c(%d)", alignment->dp_trace[i][j], alignment->dp[i][j]);
#endif
            if (i > 0 && i < end) {
                len++;
                edits++;
		if (read[j-1] == 'C')
		  insC++;
            }
            j--;
        } else {
            printf("Undefined dp_trace %d %d\n", i, j);
            alignment->ok = 0;
            return end;
        }
    }

#ifdef DEBUG_ALIGN
    printf("\nTracing done. %d %d %d %d\n", edits, len, insC, delC);
#endif

    int totC = 0;
    for (j = 0; j < m; j++) {
      if (read[j] == 'C')
	totC++;
    }
      
    
    if (edits > len * max_error_rate || delC > 6 || insC >= 0.6*totC) {
#ifdef DEBUG_ALIGN
      printf("Too many errors in read %d %d %d %d %d\n", pos[p].read, edits, (int) (len * max_error_rate), insC, delC);
#endif
      pos[p].aligned = 0;
      return end;
    }


    /* Traceback in the dp array */
    i = end;
    j = m;
    len = 0;
    edits = 0;
    old_end = end;

    while (i >= 0 && alignment->dp_trace[i][j] == 'J') {
#ifdef DEBUG_ALIGN
        printf("%c", alignment->dp_trace[i][j]);
#endif
        i--;
    }

    while (j > 0) {
        if (alignment->dp_trace[i][j] == 'D' ||
                alignment->dp_trace[i][j] == 'M') {
#ifdef DEBUG_ALIGN
            printf("%c", alignment->dp_trace[i][j]);
#endif
            er[len] = read[j - 1];
            if (read[j - 1] == 0) {
                printf("Whooot! %d %d\n", j - 1, m);
            }
            len++;
            if (alignment->dp_trace[i][j] == 'D')
                edits++;
            i--;
            j--;
        } else if (alignment->dp_trace[i][j] == 'J') {
#ifdef DEBUG_ALIGN
            printf("%c", alignment->dp_trace[i][j]);
#endif
            er[len] = '-';
            len++;
            if (j > 0 && j < m && alignment->consensus[i - 1] != '-') {
                edits++;
            }
            i--;
        } else if (alignment->dp_trace[i][j] == 'I') {
#ifdef DEBUG_ALIGN
            printf("%c", alignment->dp_trace[i][j]);
#endif
            er[len] = read[j - 1];
            if (read[j - 1] == 0) {
                printf("Whooot? %d %d\n", j - 1, m);
            }
            len++;
            if (i > 0 && i < old_end) {
                edits++;
            }
            if (i == 0) {
                alignment->offset--;
            }
            j--;

            /* Insertion in the read. Make space for the insertion in consensus*/
            for (k = end; k >= i; k--) {
                alignment->contigA[k + 1] = alignment->contigA[k];
                alignment->contigC[k + 1] = alignment->contigC[k];
                alignment->contigG[k + 1] = alignment->contigG[k];
                alignment->contigT[k + 1] = alignment->contigT[k];
                alignment->contigN[k + 1] = alignment->contigN[k];
            }
            alignment->contigA[i] = 0;
            alignment->contigC[i] = 0;
            alignment->contigG[i] = 0;
            alignment->contigT[i] = 0;
            alignment->contigN[i] = 0;
            end++;
            /* add the insertion to all aligned reads */
            for (k = 0; k < p; k++) {
                if (pos[k].aligned) {
                    if (pos[k].pos >= i) {
                        pos[k].pos++;
                    } else if (pos[k].pos < i && pos[k].pos + (signed) strlen(pos[k].edited_read) > i) {
                        alignment->contigN[i]++;
                        m2 = strlen(pos[k].edited_read);
			for (kk = m2 + 1; kk >= 0; kk--) {
			  if (pos[k].pos + kk > i) {
			    pos[k].edited_read[kk] =
			      pos[k].edited_read[kk - 1];
			  } else if (pos[k].pos + kk == i) {
			    pos[k].edited_read[kk] = '-';
			  }
                        }
                    }
                }
            }
        } else {
            printf("Undefined dp_trace %d %d\n", i, j);
            alignment->ok = 0;
            return end;
        }
    }

#ifdef DEBUG_ALIGN
    printf("\nTracing done.\n");
#endif

    pos[p].pos = i;
    i--;

#ifdef DEBUG_ALIGN
    printf("Read position: %d -> %d\n", i, i + len);
    er[len] = '\0';
    printf("ER: %s\n", er);
#endif

    /* Update alignment counts */
    for (j = 0; j < len; j++) {
        switch (er[j]) {
            case 'A':
            case 'a':
                alignment->contigA[i + len - j]++;
                break;
            case 'C':
            case 'c':
                alignment->contigC[i + len - j]++;
                break;
            case 'G':
            case 'g':
                alignment->contigG[i + len - j]++;
                break;
            case 'T':
            case 't':
                alignment->contigT[i + len - j]++;
                break;
            case 'N':
            case 'n':
            case '-':
                alignment->contigN[i + len - j]++;
                break;
            default:
                printf("Trash in er %d %d/%d\n", er[j], j, len);
        }
    }

    for (j = 0; j < len; j++) {
      pos[p].edited_read[j] = er[len - j - 1];
    }
    pos[p].edited_read[len] = '\0';

    pos[p].aligned = 1;
    return end;
}

/**
 * Compute a multiple alignment.
 * - alignment: The alignment data structure.
 * - reads: Array of reads
 * - subset: input: read ids, and approximate positions for alignment (not used at the moment)
             output: read ids, positions, and read with inserted gaps
 * - size: the number of reads to align
 * Returns the number of reads aligned. This may be smaller than size if reads are specified
 * multiple times in the input data.
 */
int multi_align(align *alignment, char **reads, read_pos *subset, int size,
        double max_error_rate, int max_aligned_reads,
        int match_reward, int mm_penalty, int gap_penalty, double bad_align_constant) {
    int i, j;
    int end;

    /* Initialize alignment */
    alignment->len = 0;
    alignment->offset = 0;
    alignment->ok = 1;

    if (size > max_aligned_reads) {
#ifdef DEBUG_ALIGN
        printf("Too many reads to align %d\n", size);
#endif

        alignment->ok = 0;
        return 0;
    }

    /* Compute alphabet skewness */
    int numC = 0, numA = 0;
    int totC[size];
    for (i = 0; i < size; i++) {
      totC[i] = 0;
      for (j = 0; j < (int) strlen(reads[subset[i].read]); j++) {
	switch (reads[subset[i].read][j]) {
	case 'A':
	case 'a':
	  numA++;
	  break;
	case 'C':
	case 'c':
	  totC[i]++;
	  numC++;
	  break;
	}
      }
    }
    int alph_skew = (int) ((double) (numA + numC) / (double) numC + 0.5);
    if (alph_skew < 1)
        alph_skew = 1;

#ifdef DEBUG_ALIGN
    printf("Alphabet skew: %d (numC: %d, numA: %d)\n", alph_skew, numC, numA);
#endif
    
    /* Form the alignment */
    end = 0;
    for (i = 0; i < size && alignment->ok; i++) {
      int edits = align_two_reads(alignment, reads, 0, i, match_reward, mm_penalty, gap_penalty, alph_skew);
      if (edits > 2*bad_align_constant * MIN(totC[0], totC[i])) {
	subset[i].aligned = 0;
      } else {
	j = align_read(alignment, reads, subset, end, i, max_error_rate,
		       match_reward, mm_penalty, gap_penalty, alph_skew);
	if (j > 0)
	  end = j;
      }
      if (!alignment->ok)
	return 0;

#ifdef DEBUG_ALIGN
      printf("Consensus ends at %d\n", end);
      printf("Consensus offset is %d\n", alignment->offset);
#endif
    }

    alignment->len = end;

    update_consensus(alignment, 0, end);

#ifdef DEBUG_ALIGN
    printf("Consensus ends at %d\n", alignment->len);
#endif

    // Mark badly aligned Rmaps
    for (i = 1; i < size; i++) {
      if (subset[i].aligned) {
        int delC = 0;
        int insC = 0;
        int totC = 0;
        for (j = 0; j < (int)strlen(subset[i].edited_read); j++) {
            if (subset[i].edited_read[j] == 'C')
                totC++;
            if (alignment->consensus[subset[i].pos + j] == 'C' && subset[i].edited_read[j] != 'C') {
                insC++;
            }
            if (alignment->consensus[subset[i].pos + j] != 'C' && subset[i].edited_read[j] == 'C') {
                delC++;
            }
        }
        //printf("%d: insC: %d, delC: %d ", (int) subset[i].read, insC, delC);
        if (delC > 3 || insC >= bad_align_constant * totC) {
            subset[i].aligned = 0;
            //printf("DELETED\n");
        } else {
            //printf("\n");
        }
      }
    }


    // Compact the subset by removing unaligned reads
    i = 0;
    for (j = 0; j < size; j++) {
        if (subset[j].aligned) {
            if (i != j) {
                subset[i].read = subset[j].read;
                subset[i].pos = subset[j].pos;
                subset[i].aligned = subset[j].aligned;
                strcpy(subset[i].edited_read, subset[j].edited_read);
                subset[i].edits = subset[j].edits;
            }
            i++;
        }
    }

    return i;
}

/* Get the consensus sequence of an alignment */
int get_consensus(align *alignment, char **cons) {
    *cons = alignment->consensus;
    return alignment->len;
}
