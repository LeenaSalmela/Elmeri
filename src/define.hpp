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
#ifndef DEFINE_H
#define DEFINE_H

#include <stdint.h>

typedef unsigned char uchar;
typedef unsigned long ulong;

#define MIN(a,b) (((a)<(b))?(a):(b))

#define ALIGN_SIZE 64 /* Maximum number of Rmaps used for a multiple alignment */

#define BLOCK_SIZE 2000 /* Block size for converting Rmaps to point patterns */

#define GOOD_ALIGN 5  /* Minimum number of Rmaps that must remain to consider an alignment for correction */

#define BAD_ALIGN_CONSTANT 0.4  /* Maximum proportion of cut sites that can be added to an Rmap */

#define MAX_READ_LENGTH 4000 /* Maximum length of an Rmap */

#define MAX_CONTIG_LEN 10000 /* Maximum length of alignment */

#define NUM_THREADS 8 /* Number of threads to use */

/* Scoring scheme for alignments (gaps in the beginning are free) */
#define GAP_PENALTY 3
#define MM_PENALTY 2
#define MATCH_REWARD 2

#endif
