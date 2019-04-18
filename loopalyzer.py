#!/usr/bin/env python

################################################################################
###  LOOPalyzer : Analyzer for LOOPER polymer libraries
###  Written by Wayland Yeung (2018)
###  University of Georgia
################################################################################

from __future__ import division
import os
import sys
import argparse
import operator
import warnings
from itertools import product
from copy import deepcopy
from scipy.stats import entropy
import numpy as np
from Bio import BiopythonWarning
from Bio import SeqIO

# each analysis is an instance of this class
class Loopalyzer:
    def __init__(self, pr_l, pr_r, format, codons):
        self.pr_l = pr_l.upper()
        self.pr_r = pr_r.upper()
        self.format = format.upper()
        self.codons = codons

        self.expected_readlen = self.codons*len(self.format)
        self.formatlen = len(self.format)
        ########################################################################

        mod_ind  = [n for n, nt in enumerate(self.format) if nt=='X']
        ran_ind  = [n for n, nt in enumerate(self.format) if nt=='N']
        constant = [n for n, nt in enumerate(self.format) if nt in list('TCGA')]
        self.mod_ind, self.ran_ind, self.constant = mod_ind, ran_ind, constant

        # list all possible error messages
        self.errors = {'left_primer': 0, 'right_primer': 0, 'both_primers': 0}
        for n in [str(a+1) for a in range(self.codons)]:
            self.errors['codon'+n] = 0

        # all possible combinations for modifications
        mod_all = []
        for mod in [''.join(r) for r in product('AGCT', repeat=len(mod_ind))]:
            temp = [0]*len(self.format)
            for n,p in enumerate(mod_ind):
                temp[p] = mod[n]
            for p in ran_ind:
                temp[p] = 'N'
            for p in constant:
                temp[p] = self.format[p]
            mod_all.append(''.join(temp))
        mod_all = dict([(c, 0) for c in mod_all])

        # all possible combinations for codon
        ran_all = []
        for mod in [''.join(r) for r in product('AGCT', repeat=len(mod_ind)+len(ran_ind))]:
            temp = [0]*len(self.format)
            for n,p in enumerate(sorted(mod_ind+ran_ind)):
                temp[p] = mod[n]
            for p in constant:
                temp[p] = self.format[p]
            ran_all.append(''.join(temp))
        ran_all = dict([(c, 0) for c in ran_all])

        # one combination set for each codon position
        self.mod_all = dict([(n+1, deepcopy(mod_all)) for n in range(self.codons)])
        self.ran_all = dict([(n+1, deepcopy(ran_all)) for n in range(self.codons)])

        # print self.constant # constant nt index
        # print self.ran_ind  # random nt index
        # print self.mod_ind  # mod nt index
        # print self.errors
        # print self.mod_all
        # print self.ran_all
        self.constant_nt = self.format[self.constant[0]]

    def show_codons(self):
        for n in sorted(self.mod_all):
            print n, ':', self.mod_all[n]
        # for n in sorted(self.ran_all):
        #     print n, ':', self.ran_all[n]

    def update_codoncount(self, codons, seqcount):
        for n, c in enumerate(codons):
            r = list(c)
            for i in self.ran_ind:
                r[i] = 'N'
            r = ''.join(r)
            # print n+1, c, r, seqcount
            self.ran_all[n+1][c] += seqcount
            self.mod_all[n+1][r] += seqcount
        # self.show_codons()

    def get_aggregates(self):
        ran_agg, mod_agg = {}, {}
        for k in self.mod_all[1]:
            mod_agg[k] = 0
        for k in self.ran_all[1]:
            ran_agg[k] = 0
        for k in self.mod_all:
            for l in mod_agg:
                mod_agg[l] += self.mod_all[k][l]
        for k in self.ran_all:
            for l in ran_agg:
                ran_agg[l] += self.ran_all[k][l]
        # print mod_agg
        # print ran_agg
        return mod_agg, ran_agg

    def notate_errors(self, sequence):
        error, read = '', ''
        if self.pr_l not in sequence:
            error+='l,'
        if self.pr_r not in sequence:
            error+='r,'
        if error == '':
            read = sequence[[n for n in xrange(len(sequence)) if sequence.find(self.pr_l, n) == n][0]+len(self.pr_l):[n for n in xrange(len(sequence)) if sequence.find(self.pr_r, n) == n][-1]]
            lread = len(read)
            if lread == self.expected_readlen:
                constants = [read[n+self.constant[0]] for n in range(0, self.expected_readlen, len(self.format))]
                error+=','.join(['c'+str(n+1) for n, nt in enumerate(constants) if nt != self.constant_nt])
            else:
                error+='n'+str(lread)
            pass
        if error == '':
            error = 'p'
        else:
            error = error.rstrip(',')
        return error, read

    def enumerate_codons(self, sequence, error, seqcount):
        if error != 'p':
            return ''
        # print sequence
        codons = [sequence[i:i+self.formatlen] for i in range(0, self.expected_readlen, self.formatlen)]
        self.update_codoncount(codons, seqcount)
        return codons

    def do(self, sequence, seqcount):
        error, read = self.notate_errors(sequence)
        codons = self.enumerate_codons(read, error, seqcount)
        return error, codons

# parses the options file
def parse_options(file):
    with open(file, 'r') as r:
        dic = [l.split('#')[0] for l in r.read().split('\n')]
        dic = [[a.strip() for a in l.split('=')] for l in dic if not l.isspace() and l !='' ]
        dic = dict(dic)
    assert 'pr_l' in dic.keys()
    assert 'pr_r' in dic.keys()
    assert 'format' in dic.keys()
    assert 'codons' in dic.keys()
    dic['codons'] = int(dic['codons'])
    if 'sw' in dic.keys():
        assert dic['sw'].lower() == 'true' or dic['sw'].lower() == 'false'
        if dic['sw'].lower() == 'true':
            dic['sw'] = True
        if dic['sw'].lower() == 'false':
            dic['sw'] = False
    return dic

# generate the counts file
def count_sequence(input):
    output = input[:-6]+'.count.tsv'

    if input[-6:] != '.fasta' and input[-6:] != '.fastq':
        sys.exit('Input file requires extension: FASTA or FASTQ')
    if input[-6:] == '.fasta':
        generator = SeqIO.parse(input, "fasta")
        sys.stderr.write('READING : %s\n' % input)
    if input[-6:] == '.fastq':
        generator = SeqIO.parse(input, "fastq")
        sys.stderr.write('READING : %s\n' % input)
        sys.stderr.write('While this program accepts FASTQ, I recommend FASTA for faster processing\n')

    aggregate = {}
    for n, record in enumerate(generator):
        if n%200000 == 0 and n!=0:
            sys.stderr.write('COUNTING: %s\n' % n)
        if record.seq not in aggregate:
            aggregate[str(record.seq)]  = 1
        else:
            aggregate[str(record.seq)] += 1
    sys.stderr.write('COUNTING: %s\n' % n)

    sorted_aggregate = sorted(aggregate.items(), key=operator.itemgetter(1), reverse=True)
    total, pad = sum([n[1] for n in sorted_aggregate]), len(str(sorted_aggregate[0][1]))
    # sys.stderr.write(pad)

    with open(output, 'w') as w:
        sys.stderr.write('WRITING : %s\n\n' % output)
        for seq, n in sorted_aggregate:
            w.write('%s\t%s\t%s\n' % (str(n).rjust(pad), format(n/total, '.9f'), seq))
    return output

# jenson shannon divergence of two probability vectors
def divergence(p, q):
    p, q = np.array(p).astype(float), np.array(q).astype(float)
    p /= p.sum()
    q /= q.sum()
    m = (p + q) / 2
    return (entropy(p, m) + entropy(q, m)) / 2

# start here
def main():
    warnings.simplefilter('ignore', BiopythonWarning)

    # handle the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-op',
                        help='mandatory options file')
    # parser.add_argument('-count', action='store_true',
    #                     help='Only generate count file')
    parser.add_argument('-r', action='store_true',
                        help='Regenerate intermediate files')
    parser.add_argument('-b', default=None,
                        help='Reference codons distribution for calculating divergence')
    args, unknown = parser.parse_known_args()
    input = unknown[0]
    if len(unknown) != 1:
        sys.exit('Please provide a single fastq or fasta file as input.')

    # # if you want to only generate the counts file and nothing else
    # if args.count:
    #     count_sequence(input)
    #     sys.exit()

    # parse the options file and output them
    option = parse_options(args.op)
    sys.stderr.write('CURRENT OPTIONS:\n  '+'\n  '.join([key.ljust(9, ' ')+str(option[key]) for key in sorted(option)])+'\n\n')

    # generate the counts file unless it already exists
    counts_file = input[:-6]+'.count.tsv'
    if not os.path.isfile(counts_file) or args.r:
        count_sequence(input)
    else:
        sys.stderr.write('Detected counts file: %s\nProceeding with this file; to disable this behavior use the -r flag\n\n' % counts_file)

    # analyze the counts file
    analysis = Loopalyzer(option['pr_l'], option['pr_r'], option['format'], option['codons'])
    reads_file = input[:-6]+'.reads.tsv'
    sys.stderr.write('Writing sequence analysis to: %s\n' % reads_file)
    w_readsfile = open(reads_file, 'w')
    with open(counts_file, 'r') as r:
        for l in r:
            seqcount, seqportion, seq = l.rstrip('\n').split('\t')
            error, codons = analysis.do(seq, int(seqcount))
            w_readsfile.write('%s\t%s\t%s\n' % (seqcount, error, ' '.join(codons)))
    w_readsfile.close()

    # write summary
    codon_file = input[:-6]+'.codons.txt'
    sys.stderr.write('Writing codons analysis to  : %s\n' % codon_file)
    with open(codon_file, 'w') as w:
        # used to both tables
        modall, ranall = analysis.get_aggregates()

        # table for modifications
        modstr = [['']+['POS'+str(n) for n in sorted(analysis.mod_all.keys())]+['TOTAL']]
        for codon in sorted(analysis.mod_all[1].keys()):
            modstr.append([codon]+[str(analysis.mod_all[pos][codon]) for pos in sorted(analysis.mod_all.keys())]+[str(modall[codon])])
        pad = max([max([len(y) for y in x]) for x in modstr])
        mod_table = ''.join([' '.join([y.ljust(pad, ' ') for y in x])+'\n' for x in modstr])

        # table for codons
        modstr = [['']+['POS'+str(n) for n in sorted(analysis.ran_all.keys())]+['TOTAL']]
        for codon in sorted(analysis.ran_all[1].keys()):
            modstr.append([codon]+[str(analysis.ran_all[pos][codon]) for pos in sorted(analysis.ran_all.keys())]+[str(ranall[codon])])
        pad = max([max([len(y) for y in x]) for x in modstr])
        ran_table = ''.join([' '.join([y.ljust(pad, ' ') for y in x])+'\n' for x in modstr])

        # list of errors
        errors = {'NO ERRORS': 0, 'LEFT PRIMER': 0, 'RIGHT PRIMER': 0, 'BOTH PRIMERS': 0}
        codons_table = dict([("CODON "+str(n+1), 0) for n in range(option['codons'])])
        codon_errors = {"LEFT FIRST":deepcopy(codons_table), "RIGHT FIRST":deepcopy(codons_table), "TOTAL":deepcopy(codons_table)}
        length_errors = {}
        with open(reads_file, 'r') as r:
            for l in r:
                count, error, seq = l.split('\t')
                if error == 'p':
                    errors['NO ERRORS'] += int(count)
                elif error == 'l,r':
                    errors['BOTH PRIMERS'] += int(count)
                elif error == 'l':
                    errors['LEFT PRIMER'] += int(count)
                elif error == 'r':
                    errors['RIGHT PRIMER'] += int(count)
                elif error[0] == 'n':
                    e = error[1:]
                    if e not in length_errors:
                        length_errors[e] = 0
                    length_errors[e] += int(count)
                elif error[0] == 'c':
                    e = [int(n[1:]) for n in error.split(',')]
                    codon_errors['LEFT FIRST']['CODON '+str(min(e))] += int(count)
                    codon_errors['RIGHT FIRST']['CODON '+str(max(e))] += int(count)
                    for h in e:
                        codon_errors['TOTAL']['CODON '+str(h)] += int(count)
        w.write(("NO ERRORS             : %s\n" % errors['NO ERRORS'])  +
                ("LEFT PRIMER MISMATCH  : %s\n" % errors['LEFT PRIMER'])+
                ("RIGHT PRIMER MISMATCH : %s\n" % errors['RIGHT PRIMER'])+
                ("BOTH PRIMERS MISMATCH : %s\n" % errors['BOTH PRIMERS']))
        if length_errors == {}:
            w.write("\nLENGTH ERRORS         : 0\n")
        else:
            w.write("\nLENGTH ERRORS         \n")
            for k in sorted(length_errors.keys()):
                w.write("LENGTH %s: %s\n" % (str(k).ljust(15,' '), length_errors[k]))

        keys = ["CODON "+str(n+1) for n in range(option['codons'])]
        w.write('\nPOLYMERIZATION ERRORS (TOTAL)\n')
        for k in keys:
            w.write('  %s : %s\n' % (k, codon_errors['TOTAL'][k]))
        w.write('\n')
        w.write('POLYMERIZATION ERRORS (FROM LEFT)\n')
        for k in keys:
            w.write('  %s : %s\n' % (k, codon_errors['LEFT FIRST'][k]))
        w.write('\n')
        w.write('POLYMERIZATION ERRORS (FROM RIGHT)\n')
        for k in keys:
            w.write('  %s : %s\n' % (k, codon_errors['RIGHT FIRST'][k]))
        w.write('\n')

        # write the tables to file
        w.write('###########################################\n### TABLE OF MODIFICATIONS PER POSITION ###\n###########################################\n')
        w.write(mod_table)
        w.write('\n####################################\n### TABLE OF CODONS PER POSITION ###\n####################################\n')
        w.write(ran_table)

        # divergence calculations
        if args.b == None:
            sys.stderr.write('\nNo reference distribution provided, skipping divergence calculation\n')
        elif not os.path.isfile(args.b):
            sys.stderr.write('\nReference distribution file not found: %s\n' % args.b)
        else:
            with open(args.b, 'r') as r:
                ref = [l.split() for l in r.read().split('\n') if len(l.split()) == 2]
                ref = dict([[c, float(n)] for c, n in ref])
            if set([e in ref for e in set(modall.keys())]) != set([True]):
                sys.stderr.write('\nReference distribution file is missing some codons\n' % args.b)
            else:
                divergence_file = input[:-6]+'.divergence.txt'
                sys.stderr.write('Writing codons divergence to: %s\n' % divergence_file)
                with open(divergence_file, 'w') as w:
                    ck, divg = sorted(modall.keys()), {}
                    reference_distribution = [ref[c] for c in ck]
                    divg['OVERALL'] = divergence(reference_distribution, [modall[c] for c in ck])
                    for n in analysis.mod_all:
                        divg['CODON '+str(n)] = divergence(reference_distribution, [analysis.mod_all[n][c] for c in ck])
                    for k in sorted(divg):
                        w.write('%s :\t%0.5f\n' % (k, divg[k]))

if __name__ == "__main__":
    main()








#
