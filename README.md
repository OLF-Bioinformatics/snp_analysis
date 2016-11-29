# snp_analysis

This is a rewrite of the SNP pipeline from https://github.com/USDA-VS/snp_analysis.

WARNING: Proper validation has not been done yet!

## Description
TODO

## Dependencies
TODO

## Usage

1. Set and run email_loopFiles2.sh
2. Set and run vcftofasta4.sh

- `vcftofasta2.sh` should give identical results than the original script. Most of the code has been made more portable.
- `vcftofasta3.sh` an hybrid between version 2 and 4. Not fully tested.
- `vcftofasta4.sh` most of the working code has been ported to Perl and Python. Huge performance improvement over the original script. Gives same results (SNPs) as original script. Previous ver
