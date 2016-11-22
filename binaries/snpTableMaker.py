#!/usr/bin/env python
import os
import random


__author__ = 'duceppemo'


class SnpTableMaker(object):
    """
    Everything is ran inside the class because data structures have to be
    shared across parent and child process during multi threading
    """

    def __init__(self, args):
        import os
        import sys
        import glob
        import multiprocessing

        # Define variables based on supplied arguments
        self.args = args
        self.ref = args.ref
        if not os.path.isfile(self.ref):
            sys.exit('Supplied reference genome file does not exists.')
        self.vcf = args.vcf
        if not os.path.isdir(self.vcf):
            sys.exit('Supplied VCF folder does not exists.')
        self.minQUAL = args.minQUAL
        if not isinstance(self.minQUAL, (int, long)):
            sys.exit('minQual value must be an integer')
        self.ac1_report = args.ac1
        self.section4 = args.section4
        self.output = args.output
        if not os.path.isdir(self.output):
            os.makedirs(self.output)
        self.table = args.table

        # number of threads to use = number of cpu
        self.cpus = int(multiprocessing.cpu_count())

        # create dictionaries to hold data
        self.refgenome = dict()
        self.vcfs = dict()
        self.ac1s = dict()
        self.ac2s = dict()
        self.allac2 = dict()
        self.finalac1 = dict()
        self.fastas = dict()
        self.counts = dict()
        self.informative_pos = dict()

        # create a list of vcf files in vcfFolder
        self.vcfList = list()
        for filename in glob.glob(os.path.join(self.vcf, '*.vcf')):
            self.vcfList.append(filename)

        # run the script
        self.snp_table_maker()

    def snp_table_maker(self):
        self.parse_ref()
        self.parse_vcf()
        self.find_ac1_in_ac2()
        self.write_ac1_report()
        self.get_allele_values()
        self.get_informative_snps()
        self.count_snps()
        self.write_fasta()
        self.write_root()
        self.write_snp_table()

    def parse_ref(self):
        from Bio import SeqIO

        print '  Parsing reference genome'

        fh = open(self.ref, "rU")
        self.refgenome = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
        fh.close()

    def parse_vcf(self):
        import sys

        print '  Parsing VCF files'

        for samplefile in self.vcfList:
            sample = os.path.basename(samplefile).split('.')[0]  # get what's before the first dot
            self.vcfs[sample] = dict()
            with open(samplefile, 'r') as f:  # open file
                for line in f:  # read file line by line
                            line = line.rstrip()  # chomp -> remove trailing whitespace characters
                            if line:  # skip blank lines or lines with only whitespaces
                                if line.startswith('##'):  # skip comment lines
                                    continue
                                elif line.startswith('#CHROM'):
                                    sample_name = line.split("\t")[9]
                                    if sample_name != sample:
                                        sys.exit('File name and sample name inside VCF file are different: %s'
                                                 % samplefile)
                                else:
                                    # chrom, pos, alt, qual = [line.split()[i] for i in (0, 1, 4, 5)]
                                    chrom = line.split()[0]
                                    pos = int(line.split()[1])
                                    alt = line.split()[4]
                                    qual = line.split()[5]  # string -> needs to be converted to integer
                                    if qual != '.':
                                        try:
                                            qual = float(qual)
                                        except ValueError:
                                            qual = int(qual)
                                    else:
                                        continue  # skip line
                                    ac = line.split()[7].split(';')[0]

                                    # http://www.saltycrane.com/blog/2010/02/python-setdefault-example/
                                    self.vcfs.setdefault(sample, {}).setdefault(chrom, {}).setdefault(pos, [])\
                                        .append(alt)

                                    if ac == 'AC=1' and qual > self.args.minQUAL:
                                        self.ac1s.setdefault(sample, {}).setdefault(chrom, []).append(pos)

                                    elif ac == 'AC=2' and qual > self.args.minQUAL:
                                        self.ac2s.setdefault(sample, {}).setdefault(chrom, []).append(pos)

                                        # This is equivalent, but faster?
                                        try:
                                            if pos not in self.allac2[chrom]:  # only add is not already present
                                                self.allac2.setdefault(chrom, []).append(pos)
                                        except KeyError:  # chromosome does not exist in dictionary
                                            self.allac2.setdefault(chrom, []).append(pos)

                                        # This works
                                        # if chrom in self.allac2:
                                        #     if pos in self.allac2[chrom]:
                                        #         pass
                                        #     else:
                                        #         self.allac2.setdefault(chrom, []).append(pos)
                                        # else:
                                        #     self.allac2.setdefault(chrom, [])

    def find_ac1_in_ac2(self):
        print '  Finding AC=1/AC=2 positions'

        if isinstance(self.ac1s, dict):  # check if it's a dict before using .iteritems()
            for sample, chromosomes in self.ac1s.iteritems():
                if isinstance(chromosomes, dict):  # check for dict
                    for chrom, positions in chromosomes.iteritems():
                        if isinstance(positions, list):  # check for list
                            for pos in positions:
                                if pos in self.allac2[chrom]:  # check ac1 in ac2
                                    self.finalac1.setdefault(sample, {}).setdefault(chrom, []).append(pos)

    def write_ac1_report(self):
        print "  Writing AC=1/AC=2 report to file"

        # free up resources not needed anymore
        self.ac1s.clear()

        fh = open(self.ac1_report, 'w')
        if isinstance(self.finalac1, dict):
            for sample, chromosomes in sorted(self.finalac1.iteritems()):
                if isinstance(chromosomes, dict):
                    for chrom, positions in sorted(chromosomes.iteritems()):
                        if isinstance(positions, list):
                            fh.write("{}\nAC=1 is also found in AC=2 in chromosome {}".format(sample, chrom) +
                                     " at position(s): " + ', '.join(map(str, positions)) + "\n\n")
        fh.close()

    def get_allele_values(self):
        print '  Getting allele values'

        for sample in self.ac2s:
            for chrom in self.ac2s[sample]:
                for pos in self.allac2[chrom]:
                    # if in AC=2 for that sample
                    if pos in self.ac2s[sample][chrom]:
                        allele = ''.join(self.vcfs[sample][chrom][pos])  # convert list to string
                    else:
                        try:  # use a try here because some samples are not in finalac1
                            # if in AC=1 for that sample, but also in AC=2 in other sample
                            if pos in self.finalac1[sample][chrom]:
                                allele = ''.join(self.vcfs[sample][chrom][pos])  # convert list to string
                            else:
                                allele = self.refgenome[chrom].seq[pos - 1]
                        except KeyError:
                            allele = self.refgenome[chrom].seq[pos - 1]

                    self.fastas.setdefault(sample, {}).setdefault(chrom, {}).setdefault(pos, []).append(allele)

                    # Track all alleles for each position
                    try:
                        if allele not in self.counts[chrom][pos]:
                            self.counts.setdefault(chrom, {}).setdefault(pos, []).append(allele)
                    except KeyError:
                        self.counts.setdefault(chrom, {}).setdefault(pos, []).append(allele)

    def get_informative_snps(self):
        """SNPs position that have at least one different ALT allele within all the samples"""

        print '  Getting informative SNPs'

        # free up resources not needed anymore
        self.ac2s.clear()
        self.allac2.clear()
        self.finalac1.clear()

        # need to get the positions in the same order for all the sample (sort chrom and pos)
        for sample in self.fastas:
            for chrom in sorted(self.fastas[sample]):
                for pos in sorted(self.fastas[sample][chrom]):
                    if len(self.counts[chrom][pos]) > 1:  # if more that one ALT allele, keep it
                        allele = ''.join(self.fastas[sample][chrom][pos])  # convert list to string

                        # check if allele is empty
                        if allele:
                            self.informative_pos.setdefault(sample, {}).setdefault(chrom, {})\
                                .setdefault(pos, []).append(''.join(allele))
                        else:
                            print "No allele infor for {}, {}:{}".format(sample, chrom, pos)

    def count_snps(self):
        print '  Counting SNPs'

        # free up resources not needed anymore
        self.counts.clear()

        # All samples should have the same number of informative SNPs
        # so any can be used to get the stats
        randomsample = random.choice(self.informative_pos.keys())

        filteredcount = 0
        informativecount = 0

        # Account for multiple chromosome
        for chrom in self.fastas[randomsample]:
            filteredcount += len(self.fastas[randomsample][chrom])  # number of positions
            informativecount += len(self.informative_pos[randomsample][chrom])

        # print to screen
        print "\nTotal filtered SNPs: {}".format(filteredcount)
        print "Total informative SNPs: {}\n".format(informativecount)

        # write to file
        fh = open(self.section4, "a")  # append mode
        fh.write("Total filtered SNPs: {}\n".format(filteredcount))
        fh.write("Total informative SNPs: {}\n\n".format(informativecount))
        fh.close()

    def write_fasta(self):
        print '  Writing sample fasta files'

        # free up resources not needed anymore
        self.fastas.clear()

        # Create output folder for fasta files
        if not os.path.exists(self.output):
            os.makedirs(self.output)

        if isinstance(self.informative_pos, dict):
            for sample, chromosomes in sorted(self.informative_pos.iteritems()):
                samplepath = os.path.join(self.output, sample + '.fas')
                fh = open(samplepath, 'w')
                fh.write(">{}\n".format(sample))
                if isinstance(chromosomes, dict):
                    for chrom, positions in sorted(chromosomes.iteritems()):
                        if isinstance(positions, dict):
                            for pos, allele in sorted(positions.iteritems()):
                                if isinstance(allele, list):
                                    fh.write(''.join(allele))  # convert list to text
                fh.write("\n")

    def write_root(self):
        print '  Writing root fasta file'

        rootpath = os.path.join(self.output, 'root.fas')
        randomsample = random.choice(self.informative_pos.keys())
        rootseq = list()

        fh = open(rootpath, 'w')
        if isinstance(self.informative_pos, dict):
            for chrom in self.informative_pos[randomsample]:
                for pos in sorted(self.informative_pos[randomsample][chrom]):
                    rootseq.append(self.refgenome[chrom].seq[pos - 1])
        fh.write(">root\n" + "{}\n".format(''.join(rootseq)))

    def write_snp_table(self):
        print '  Writing SNP table'

        fh = open(self.table, 'w')

        randomsample = random.choice(self.informative_pos.keys())

        ref_pos = list()
        ref_call = list()

        # reference
        if isinstance(self.informative_pos, dict):
            for chrom in self.informative_pos[randomsample]:
                for pos in sorted(self.informative_pos[randomsample][chrom]):
                    ref_pos.append(''.join(chrom) + '-' + str(pos))
                    ref_call.append(self.refgenome[chrom].seq[pos - 1])

        fh.write("reference_pos\t{}\n".format("\t".join(ref_pos)))
        fh.write("reference_call\t{}\n".format("\t".join(ref_call)))

        # sample
        if isinstance(self.informative_pos, dict):
            for sample, chromosomes in self.informative_pos.iteritems():
                fh.write("{}".format(sample))
                if isinstance(chromosomes, dict):
                    for chrom, positions in sorted(chromosomes.iteritems()):
                        if isinstance(positions, dict):
                            for pos, allele in sorted(positions.iteritems()):
                                if isinstance(allele, list):
                                    allele = ''.join(allele)  # convert list to text
                                    fh.write("\t{}".format(allele))
                fh.write("\n")
        fh.close()


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Generate SNP table and aligned fasta files from VCF files')
    parser.add_argument('-r', '--ref', metavar='ref.fasta',
                        required=True,
                        help='reference genome used in the VCF files')

    parser.add_argument('-v', '--vcf', metavar='vcfFolder',
                        required=True,
                        help='location of the VCF files')

    parser.add_argument('-q', '--minQUAL', metavar='minQUAL', type=int,
                        required=True,
                        help='minimum QUAL value in VCF file')

    parser.add_argument('-ac1', '--ac1', metavar='AC1Report.txt',
                        required=True,
                        help='output file where positions having both AC=1 and AC=2 are reported')

    parser.add_argument('-s4', '--section4', metavar='section4.txt',
                        required=True,
                        help='output file where total filtered SNP positions and total informative SNPs are reported')

    parser.add_argument('-o', '--output', metavar='fastaOutFolder',
                        required=True,
                        help='folder where the output fasta files will be output')

    parser.add_argument('-t', '--table', metavar='fastaTable.tsv',
                        required=True,
                        help='the SNP table')

    # Get the arguments into an object
    arguments = parser.parse_args()

    SnpTableMaker(arguments)
