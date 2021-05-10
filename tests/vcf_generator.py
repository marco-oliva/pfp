#!/usr/bin/env python
# coding: utf-8


import gzip
import random
import math
import operator
from tqdm import tqdm
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import argparse
from datetime import datetime



# Parameters
class Parameters:
    r = 2
    s = 251
    n = 10000
    p = 20
    N = 0
    P = 70
    offset = 17000000
    R = '/blue/boucher/marco.oliva/tmp/generated_vcf/22.fa.gz'
    out_file = 'generated.vcf.gz'


Description = '''
VCF Generator
'''

# ### Spec
# Generate a VCF file so that:
#
# - **R** refernce path
# - **r** groups, wihithin a group samples are similar (simulating region similarities)
# - **s** samples per group, deepness
# - **n** random variation sites for each group
# - **N** variations common to all regions and samples
# - **p** out of 100 samples in a group will have a certain SNP
# - **P** same as p for N
# - **offset** of the first possible position for a variation

def main():

    parameters = Parameters()

    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-R', help='reference path', type=str, dest='reference', required=True)
    parser.add_argument('-r', help='input B file name', type=int, default=parameters.r, dest="groups")
    parser.add_argument('-s', help='samples per group', type=int, default=parameters.s, dest="samples")
    parser.add_argument('-n', help='random variation sites for each group', type=int, default=parameters.n, dest="variations")
    parser.add_argument('-N', help='variations common to all regions and samples', type=int, default=parameters.N, dest="common")
    parser.add_argument('-p', help='out of 100 samples in a group will have a certain SNP', type=int, default=parameters.p, dest="snp_prob")
    parser.add_argument('-P', help='same as p for N', type=int, default=parameters.P, dest="common_snp_probs")
    parser.add_argument('-o', help='output file', type=str, default=parameters.out_file, dest="out_file")
    parser.add_argument('--offset', help='the first possible position for a variation', type=int, default=parameters.offset, dest="offset")
    args = parser.parse_args()


    # Set parameters
    parameters.r = args.groups
    parameters.s = args.samples
    parameters.n = args.variations
    parameters.p = args.snp_prob
    parameters.N = args.common
    parameters.P = args.common_snp_probs
    parameters.offset = args.offset
    parameters.R = args.reference
    parameters.out_file = args.out_file

    # Read the reference
    reference = ""
    reference_length = 0
    with gzip.open(parameters.R, "rt") as handle:
        record = SeqIO.read(handle, "fasta")
        reference = record.seq
        reference_length = len(reference)

    # Variation Class
    class Variation:
        chromosome = 'GEN'
        pos = 0
        variation_id = ''
        ref = ''
        alt = ''
        qual = '100'
        filter_s = 'PASS'
        info = '.'
        format_s = 'GT'
        samples = list()


    # Variations list (vcf)
    variations = list()

    # Seeding random generator with time, something more random?
    random.seed(datetime.now())

    # Create samples per region map
    samples_per_region = dict()
    for i in range(parameters.r):
        samples_list = list()
        for sample_n in range(parameters.s):
            samples_list.append('region_{}_sample_{}'.format(i,sample_n))
        samples_per_region[i] = samples_list

    # For each region
    print('Generating random variations for each region')
    rand_variations_pos = list()
    for region, samples_list in samples_per_region.items():
        # Create variations
        for i in tqdm(range(parameters.n)):

            # Get variation random pos
            while True:
                pos = math.floor(random.uniform(parameters.offset, reference_length))
                if (pos not in rand_variations_pos):
                    break
            rand_variations_pos.append(pos)

            # Samples subset
            samples_sublist = random.sample(samples_list, k=(math.floor((parameters.p/100)*len(samples_list))))

            # Build variation object
            var = Variation()
            var.pos = pos
            var.ref = reference[pos]

            nucleotide_set = {'A','C','G','T','N'}
            nucleotide_set.remove(reference[pos])
            var.alt = random.sample(nucleotide_set, k=1)[0]

            var.reference_length = 1
            var.samples = samples_sublist
            var.variation_id = 'region_{}_variation_{}'.format(region, i)

            # Append variation to global variations list
            variations.append(var)

    # Add variations common to all samples
    all_samples = list()
    for region, samples_list in samples_per_region.items():
        all_samples.extend(samples_list)

    print('Building vcf structure')
    for v in range(parameters.N):
        # Build variation object
        while True:
            pos = math.floor(random.uniform(parameters.offset, reference_length))
            if (pos not in rand_variations_pos):
                break
        rand_variations_pos.append(pos)

        samples_sublist = random.sample(all_samples, k=(math.floor((parameters.P/100)*len(all_samples))))

        var = Variation()
        var.pos = pos
        var.ref = reference[pos]

        nucleotide_set = {'A','C','G','T','N'}
        nucleotide_set.remove(reference[pos])
        var.alt = random.sample(nucleotide_set, k=1)[0]

        var.reference_length = 1
        var.samples = samples_sublist
        var.variation_id = 'all_samples_variation_{}'.format(region, i)

        # Append variation to global variations list
        variations.append(var)

    # Sort Variations by position
    variations.sort(key=operator.attrgetter('pos'))

    # Variations statistics
    print('==============================')
    print('Statistics')

    positions = np.array([v.pos for v in variations])
    print('Min pos: {}'.format(np.min(positions)))
    print('Max pos: {}'.format(np.max(positions)))

    distances = np.zeros(positions.size)
    for idx, e in enumerate(positions):
        if (idx == distances.size - 1):
            break
        distances[idx] = positions[idx + 1] - positions[idx]
    distances[distances.size - 1] = distances[distances.size - 2]

    print('Dist < 40: {}'.format(len([np.where( distances < 40 )])))
    print('Min dist: {}'.format(np.min(distances)))
    print('Max dist: {}'.format(np.max(distances)))
    print('Average dist: {}'.format(np.average(distances)))
    print('==============================')

    data = distances

    # log-scaled bins
    bins = np.logspace(0, 5, 50)
    widths = (bins[1:] - bins[:-1])

    # Calculate histogram
    hist = np.histogram(data, bins=bins)
    # normalize by bin width
    hist_norm = hist[0]/widths

    # plot it!
    plt.bar(bins[:-1], hist_norm, widths)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Distance between two variations in nucleotides')
    plt.ylabel('Abundance')
    plt.savefig(parameters.out_file + '.plot.png')

    # Output VCF file
    print('Writing to file')
    with gzip.open(parameters.out_file, 'wb') as f:
        f.write('##fileformat=VCFv4.1\n'.encode())
        f.write('##reference={}\n'.format(parameters.R).encode())
        f.write('##r:\t{}\n'.format(parameters.r).encode())
        f.write('##s:\t{}\n'.format(parameters.s).encode())
        f.write('##n:\t{}\n'.format(parameters.n).encode())
        f.write('##p:\t{}\n'.format(parameters.p).encode())
        f.write('##N:\t{}\n'.format(parameters.N).encode())
        f.write('##P:\t{}\n'.format(parameters.P).encode())
        f.write('##offset:\t{}\n'.format(parameters.offset).encode())
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'.encode())

        header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format('\t'.join(all_samples))
        f.write(header.encode())

        # Write Variations
        for variation in tqdm(variations):
            variation_line = list()
            variation_line.append(variation.chromosome)
            variation_line.append(str(variation.pos))
            variation_line.append(variation.variation_id)
            variation_line.append(variation.ref)
            variation_line.append(variation.alt)
            variation_line.append(variation.qual)
            variation_line.append(variation.filter_s)
            variation_line.append(variation.info)
            variation_line.append(variation.format_s)

            # Samples matrix
            for sample in all_samples:
                if (sample in variation.samples):
                    variation_line.append('1|1')
                else:
                    variation_line.append('0|0')

            # To string
            variation_string = '\t'.join(variation_line)
            f.write(variation_string.encode())
            f.write('\n'.encode())



if __name__ == '__main__':
    main()