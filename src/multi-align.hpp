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

#ifndef MULTI_ALIGN_H
#define MULTI_ALIGN_H

#include "define.hpp"

/* Information about the reads in the alignment */
typedef struct {
    uint64_t read;      /* The id of the read. (input) */
    int pos;            /* The position of the read in the alignment.
                           (input and output, input should give an
                           estimate of the positioning of the read) */
    char aligned;       /* Set when the read has been aligned */

    char *edited_read;  /* The read with gaps added according to the alignment (output) */
    int edits;          /* Number of edit operations as compared to consensus */
} read_pos;

typedef struct {
    /* Arrays for dynamic programming */
    int dp[MAX_CONTIG_LEN+10][MAX_CONTIG_LEN+10];
    char dp_trace[MAX_CONTIG_LEN+10][MAX_CONTIG_LEN+10];

    /* Arrays for counts of bases in each column of the alignment */
    int contigA[MAX_CONTIG_LEN];
    int contigC[MAX_CONTIG_LEN];
    int contigG[MAX_CONTIG_LEN];
    int contigT[MAX_CONTIG_LEN];
    int contigN[MAX_CONTIG_LEN]; /* Count of gaps */

    char consensus[MAX_CONTIG_LEN]; /* The consensus of the alignment */
    int len;                        /* The length of consensus */

    int offset;  /* Offset of the begining of consensus with regard to read pos estimates */
    int ok;
} align;

/**
 * Compute a multiple alignment.
 * - alignment: The alignment data structure.
 * - reads: Array of reads
 * - subset: input: read ids, orientations and approximate positions for alignment
             output: read ids, orientations, positions, and read with inserted gaps
 * - size: the number of reads to align
 * Returns the number of reads aligned. This may be smaller than size if reads are specified
 * multiple times in the input data.
 */
int multi_align(align *alignment, char **reads, read_pos *subset, int size,
                double max_error_rate, int max_aligned_reads,
                int match_reward, int mm_penalty, int gap_penalty, double bad_align_constant);


/* Get the consensus sequence of an alignment */
int get_consensus(align *alignment, char **cons);

#endif
