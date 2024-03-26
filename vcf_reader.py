#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Created By  : Ardalan Naseri
# Created Date: Mon September 21 2020
# =============================================================================
"""The module is a VCF reader to parse input VCF file."""

import gzip
import random


def eff_split(string_input, array_vals, delimeter='\t'):
    counter = 0
    start_pos = 0
    end_pos = start_pos

    while start_pos < len(string_input) - 1:
        end_pos = start_pos + 1
        while end_pos < len(string_input) and string_input[end_pos] != delimeter and start_pos != end_pos:
            end_pos += 1
        array_vals[counter] = string_input[start_pos:end_pos]
        start_pos = end_pos + 1
        counter = counter + 1


class VCFReader:

    def __init__(self, vcf_input_compressed):
        self.vcf_file = gzip.open(vcf_input_compressed)
        self.samples = []
        self.done = False
        self.vals = []
        self.genome_pos = []
        self.valid = True
        self.entries_started = False
        self.inter_vals = None
        self._line = None

    def set_samples(self):
        done = False
        while not done:
            line = self.vcf_file.readline()
            if not line:
                done = True
                self.done = True
                continue
            if '#CHROM' in line:
                self.entries_started = True
                i = 9
                _values = line.replace("\n", "").split()
                while i < len(_values):
                    self.samples.append(_values[i])
                    i += 1
                self.vals = [0] * len(self.samples)
                self.inter_vals = ['0|1'] * (len(self.samples) + 9)
                done = True

    def read_next_site(self):

        site_counter = 0
        line = self.vcf_file.readline().replace("\n", "")
        self._line = line
        self.valid = True

        if not line:
            self.done = True
            self.vcf_file.close()
            return False

        if self.entries_started:
            #eff_split(line, self.inter_vals, '\t')
            self.inter_vals = line.split()
            _pos = self.inter_vals[1]
            alt = self.inter_vals[4]
            if len(alt.split(',')) > 1:
                self.valid = False
                return True
            i = 2
            while i < len(self.inter_vals) and self.inter_vals[i] != 'GT':
                i += 1
            i += 1

            if i >= len(self.inter_vals):
                self.valid = False
                return True
            tags = self.inter_vals[7]
            if len(self.inter_vals[3]) > 1 or len(self.inter_vals[4]) > 1:
                self.valid = False
                return True
            i = 9
            site_values = ''
            j = 0
            while i < len(self.inter_vals):
                site_values = self.inter_vals[i].replace("\n", '').split("|")

                if site_values[0] == '.' or len(site_values) < 2 or (len(site_values) > 1 and site_values[1] == '.'):
                    self.valid = False
                    return True

                al1 = int(site_values[0])
                al2 = int(site_values[1])

                if al1 == al2:
                    self.vals[j] = al1
                else:
                    self.vals[j] = random.randint(0, 1)
                j = j + 1
                i += 1

            self.genome_pos.append(self.inter_vals[1])
            site_counter = site_counter + 1
            return True
