# Elmeri

Elmeri is a program to correct errors in raw optical map data (Rmaps).

## Reference

L. Salmela, K. Mukherjee, S.J. Puglisi, M.D. Muggli, and C. Boucher:
Fast and accurate correction of optical mapping data with spaced seeds

## System Requirements

Elmeri has been tested on systems running Linux on a X86_64 architecture. Elmeri requires that the boost library is available.

## Installation

Copy the repository to your own machine and run make in the src directory.

## Usage

elmeri [parameters]

| Required parameters: |  |
| --- | --- |
| -i input_file     | Input file of Rmaps to be corrected |
| -o output_file    | Output file for corrected Rmaps |

| Optional parameters: |   |
| --- | --- |
| -l el            | Size of (el,k)-mers in kbp (default: 80) |
| -k k              | Minimum size of (el,k)-mers in number of fragments (default: 5) |
| -s spacing_pattern | Spacing pattern for spaced seeds.(default: 11111111110001110110010010011101001110001010010100001010011000010111100000001100) |
| -c coverage         | Maximum number of Rmaps in an alignment. If set to > 64, will be reset to 64. (default: 64) |
| -t similarity threshold | Threshold for merging similar (el,k)-mers in the index. (default: 5) |

The default parameters work well on small to medium sized genomes when
the coverage of the Rmaps is high. For a large genome such as the
human genome, setting k to a higher value (e.g. 6) to ensure that
(el,k)-mers are unique in the genome often gives better results. The
coverage parameter should be set to reflect the coverage of the Rmap
data.

## Data format

The input Rmaps should be in the following format. Each Rmap takes
three lines in the file. First line gives the name of the Rmap. The
second line has two columns for the enzyme and the rest of the columns
give the fragment sizes. The third line is empty.

If you have data in the bnx format, the bnx2valouev.py script in the
eval folder can be used to convert the data into the correct format:

python bnx2valouev.py input.bnx > input.valouev

## Example

An example data set is given in the ecoli-sample folder. This data set
is already in the correct format for elmeri. To correct the Rmaps run
elmeri as follows:

elmeri -i ecoli-sample/ecoli-2000.valouev -o ecoli-sample/ecoli-2000-corrected.valouev

## Outputting Rmaps sharing (el,k)-mers

The elmeri-index binary can be used to output pairs of Rmaps that
share (el,k)-mers:

elmeri-index input_Rmaps el output_file spacing_pattern count_threshold merge_similarity_threshold

| Parameters: |  |
| --- | --- |
| input_Rmaps     | Input file of Rmaps to be indexed |
| el             | Size of (el,k)-mers in kbp |
| output_file     | Output file for pairs of related Rmaps |
| spacing_pattern | Spacing pattern for spaced seeds.(e.g. 11111111110001110110010010011101001110001010010100001010011000010111100000001100) |
| count_threshold | Threshold for the number of shared (el,k)-mers two Rmaps have to share to be considered related. |
| merge_similarity_threshold | Threshold for merging similar (el,k)-mers in the index. |

The input Rmaps should be in the same format as for Elmeri. The output
is provided in dot graph format where the nodes are the Rmaps (labeled
by integers giving their order in the input file) and edges represent
relatedness of the Rmaps. The weight of the edges gives the number of
(el,k)-mers the Rmaps share.
