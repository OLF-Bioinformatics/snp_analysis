# snp_analysis

This is a rewrite of the SNP pipeline from https://github.com/USDA-VS/snp_analysis.

WARNING: This is still experimental.

- `vcftofasta2.sh` should give identical results than the original script. Most of the code has been made more portable.
- `vcftofasta3.sh` an hybrid between version 2 and 4. Not fully tested
- `vcftofasta4.sh` most of the code has been ported to Perl. Huge performance improvement over the original script. Performance can still be increased by parallelizing many loops (TODO).
