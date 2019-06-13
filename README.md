# Elmeri

Elmeri is a program to correct errors in raw optical map data (Rmaps).

## Reference

L. Salmela, K. Mukherjee, S.J. Puglisi, M.D. Muggli, and C. Boucher:
Fast and accurate correction of optical mapping data with spaced seeds

## System Requirements

Elmeri has been tested on systems running Linux on a X86_64 architecture.

## Installation

Copy the repository to your own machine and run make in the src directory.

## Usage

elmeri [parameters]

| Required parameters: |  |
| --- | --- |
| -i <input file>     | Input file of Rmaps to be corrected |
| -o <output file>    | Output file for corrected Rmaps |

| Optional parameters: |   |
| --- | --- |
| -l <ell>            | Size of (el,k)-mers in kbp (default: 80) |
| -k <k>              | Minimum size of (el,k)-mers in number of fragments (default: 5) |
| -s <spacing pattern> | Spacing pattern for spaced seeds.(default: 11111111110001110110010010011101001110001010010100001010011000010111100000001100) |
| -c <coverage>         | Maximum number of Rmaps in an alignment. If set to > 64, will be reset to 64. (default: 64) |
| -t <similarity threshold> | Threshold for merging similar (el,k)-mers in the index. (default: 5) |

## Example

elmeri -i input.valouev -o output.valouev

The input Rmaps should be in the following format. Each Rmap takes
three lines in the file. First line gives the name of the Rmap. The
second line has two columns for the enzyme and the rest of the columns
give the fragment sizes. The third line is empty.
