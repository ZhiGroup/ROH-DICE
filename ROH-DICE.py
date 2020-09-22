#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Created By  : Ardalan Naseri
# Created Date: Mon September 21 2020
# =============================================================================
"""The script implements cPBWT algorithms for finding ROH diplotypes using PBWT"""

import argparse
import bisect
import itertools
import random
import sys

from vcf_reader import VCFReader


def quick_select(_list, k):
    def select(list_of_elements, left, right, index):

        if right == left:
            return list_of_elements[left]

        pivot_index = random.randint(left, right)

        list_of_elements[left], list_of_elements[pivot_index] = list_of_elements[pivot_index], list_of_elements[left]

        i = left
        for j in xrange(left + 1, right + 1):
            if list_of_elements[j] < list_of_elements[left]:
                i = i + 1
                list_of_elements[i], list_of_elements[j] = list_of_elements[j], list_of_elements[i]

        list_of_elements[i], list_of_elements[left] = list_of_elements[left], list_of_elements[i]

        if index == i:
            return list_of_elements[i]
        elif index < i:
            return select(list_of_elements, left, i - 1, index)
        else:
            return select(list_of_elements, i + 1, right, index)

    if _list is None or len(_list) < 1:
        return None

    return select(_list, 0, len(_list) - 1, k)


class WeightedRandomGen(object):
    def __init__(self, weights):
        self.totals = []
        sum_t = 0

        for w in weights:
            sum_t += w
            self.totals.append(sum_t)

    def next(self):
        rnd = random.random() * self.totals[-1]
        return bisect.bisect_right(self.totals, rnd)

    def __call__(self):
        return self.next()


def max_number_of_samples(vcf_input_compressed, c=10, l=1000, t=2, output_fp='output.txt'):
    vcf_parser = VCFReader(vcf_input_compressed)
    vcf_parser.set_samples()
    f_o = open(output_fp, 'w+')
    M = len(vcf_parser.samples)
    # initialise positional prefix array
    ppa = list(range(M))
    # initialise divergence array
    div = [0] * M
    k = 0

    while vcf_parser.read_next_site():

        if vcf_parser.valid == False:
            continue
        # print vcf_parser.vals
        V = []  # for a and b
        Q = []  # for d and e
        for _ in range(t + 1):
            V.append(list())
            Q.append(list())
        pq = []
        mab = []
        for _ in range(0, t):
            pq.append(k + 1)
            mab.append(list())
        # iterate over haplotypes in reverse prefix sorted order
        i0 = 0
        ll = -1
        # current haplotype
        haplotype = vcf_parser.vals  # X[index]

        for index, match_start in zip(ppa, div):
            ll += 1

            if match_start > k - l:

                changed = False
                for p in range(0, t):
                    for q in range(p + 1, t):
                        if len(mab[p]) > 0 and len(mab[q]) > 0:  # or k == N-1:
                            changed = True
                            break
                    if changed:
                        break

                if changed and ll - i0 >= c:
                    d_min = 0
                    for ib in range(i0, ll):
                        f_o.write(vcf_parser.samples[ppa[ib]] + "\t")
                        # _concatenation.append(vcf_parser.samples[ppa[ib]])
                        if div[ib] > d_min:
                            d_min = div[ib]
                    f_o.write(str(vcf_parser.genome_pos[k]) + "\t" + str(vcf_parser.genome_pos[d_min]) + "\n")

                i0 = ll
                mab = []
                for _ in range(0, t):
                    mab.append(list())

            # allele for current haplotype
            allele = haplotype[index]

            # update intermediates
            for p in range(0, t):
                if match_start > pq[p]:
                    pq[p] = match_start

            # update intermediates
            V[allele].append(index)

            for p in range(0, t):
                if allele == p:
                    Q[allele].append(pq[allele])
                    pq[allele] = 0
                    mab[p].append(index)
                    # dd[p].append(ll)
        changed = False
        for p in range(0, t):
            for q in range(p + 1, t):
                if len(mab[p]) > 0 and len(mab[q]) > 0:  # or k== N-1:
                    changed = True
                    break
            if (changed):
                break

        if changed and M - 1 - i0 >= c:

            # _concatenation = []
            d_min = 0
            # _concatenation.append(vcf_parser.samples[ppa[i0]])
            for ib in range(i0, M - 1):
                f_o.write(vcf_parser.samples[ppa[ib]] + "\t")
                # _concatenation.append(vcf_parser.samples[ppa[ib]])
                # for ib in range(ia+1,M-1):
                if div[ib] > d_min:
                    d_min = div[ib]
            f_o.write(str(vcf_parser.genome_pos[k]) + "\t" + str(vcf_parser.genome_pos[d_min]) + "\n")

        i0 = ll
        mab = []
        for _ in range(0, t):
            mab.append(list())

        ppa = list(itertools.chain(*V))
        div = list(itertools.chain(*Q))
        k = k + 1


def max_number_of_sites(vcf_input_compressed, c=10, l=1000, t=2, output_fp='output.txt'):

    vcf_parser = VCFReader(vcf_input_compressed)
    vcf_parser.set_samples()

    M = len(vcf_parser.samples)
    # initialise positional prefix array
    ppa = list(range(M))
    # initialise divergence array
    div = [0] * M
    k = 0
    f_o = open(output_fp, 'w+')

    while vcf_parser.read_next_site():

        if not vcf_parser.valid:
            continue
        V = []  # for a and b
        Q = []  # for d and e
        for _ in range(t + 1):
            V.append(list())
            Q.append(list())
        pq = []
        mab = []
        for _ in range(0, t):
            pq.append(k + 1)
            mab.append(list())
        # iterate over haplotypes in reverse prefix sorted order
        i0 = 0
        ll = -1
        # current haplotype
        haplotype = vcf_parser.vals  # X[index]

        i0 = 0
        ll = -1
        for index, match_start in zip(ppa, div):
            # print index
            ll += 1
            if match_start > k - l:
                changed = False
                for p in range(0, t):
                    for q in range(p + 1, t):
                        if (len(mab[p]) > 0 and len(mab[q]) > 0):  # or k == N-1:
                            changed = True
                            break
                    if (changed):
                        break
                if changed and ll - i0 >= c:
                    report = True

                    if report:
                        T_d = []
                        all_d = []
                        for l in range(0, t):
                            T_d.append(list())

                        _concatenation = []
                        d_min = 0
                        al = haplotype[ppa[i0]]
                        dv_first = div[i0 + 1]

                        T_d[al].append(dv_first)
                        all_d.append(dv_first)
                        for ia in range(i0 + 1, ll):
                            d_min = 0
                            if div[ia] > d_min:
                                d_min = div[ia]
                            al = haplotype[ppa[ia]]
                            T_d[al].append(d_min)
                            all_d.append(d_min)
                        # all_d.sort()
                        skip_this = False
                        for l in range(0, t):
                            # T_d[l].sort()
                            # quick_select(T_d)
                            if len(T_d[l]) >= c and quick_select(T_d[l], c - 1) <= quick_select(all_d, c - 1):
                                skip_this = True
                                break
                        if not skip_this:
                            ia = i0
                            if (dv_first <= quick_select(all_d,
                                                         c - 1)):  # all_d [c-1]): # s_p_m[vcf_parser.samples[ppa[i0]]]
                                f_o.write(vcf_parser.samples[ppa[i0]] + "\t")

                            for ia in range(i0 + 1, ll):
                                if div[ia] <= quick_select(all_d, c - 1):  # all_d[c-1]):
                                    f_o.write(vcf_parser.samples[ppa[ia]] + "\t")
                            d_res = quick_select(all_d, c - 1)
                            f_o.write(
                                str(vcf_parser.genome_pos[k - 1]) + "\t" + str(vcf_parser.genome_pos[d_res]) + "\n")

                i0 = ll
                mab = []
                for _ in range(0, t):
                    mab.append(list())

            #            allele = haplotype[ll]
            allele = haplotype[index]

            # update intermediates
            for p in range(0, t):
                if match_start > pq[p]:
                    pq[p] = match_start

            # update intermediates
            V[allele].append(index)

            for p in range(0, t):
                if allele == p:
                    Q[allele].append(pq[allele])
                    pq[allele] = 0
                    mab[p].append(index)
                    # dd[p].append(ll)
        changed = False
        for p in range(0, t):
            for q in range(p + 1, t):
                if len(mab[p]) > 0 and len(mab[q]) > 0:  # or k== N-1:
                    changed = True
                    break
            if (changed):
                break
        if changed and M - 1 - i0 >= c:
            report = True
            for l in range(0, t):
                if mab[l] >= c:
                    report = False
                    break
            if (report):
                _concatenation = []
                d_min = 0
                f_o.write(vcf_parser.samples[ppa[i0]] + "\t")
                for ib in range(i0 + 1, M - 1):
                    f_o.write(vcf_parser.samples[ppa[ib]] + "\t")
                    # for ib in range(ia+1,M-1):
                    if div[ib] > d_min:
                        d_min = div[ib]
                f_o.write(str(vcf_parser.genome_pos[k - 1]) + "\t" + str(vcf_parser.genome_pos[d_min]) + "\n")

        i0 = ll
        mab = []
        for _ in range(0, t):
            mab.append(list())

        ppa = list(itertools.chain(*V))
        div = list(itertools.chain(*Q))
        k = k + 1


def get_options(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="ROH-DICE")
    parser.add_argument("-i", "--input", help="Input file in compressed (.gz) VCF format")
    parser.add_argument("-o", "--output", help="Your destination output file.")
    parser.add_argument("-w", "--min_samples", type=int, help="Minimum number of samples in each ROH cluster.")
    parser.add_argument("-l", "--min_length", type=int, help="Minimum number of sites in each ROH cluster.")
    parser.add_argument("-L", "--max_length", action='store_true',
                        help="Maximize the number of sites in each cluster (by default samples are maximized).")
    options = parser.parse_args(args)
    return options


options = get_options(sys.argv[1:])
max_length = False

if options.max_length:
    max_length = True

vcf_c = options.input
output_file = options.output
min_samples = options.min_samples
min_length = options.min_length

if max_length:
    max_number_of_sites(vcf_c, min_samples, min_length, 2, output_file)

else:
    max_number_of_samples(vcf_c, min_samples, min_length, 2, output_file)
