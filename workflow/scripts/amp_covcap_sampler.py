# Copyright 2022 Simon Magin.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed except according to those terms.

import random
from collections import Counter, defaultdict


class Mapping:
    def __init__(
        self,
        qname,
        qlen,
        qstart,
        qend,
        samestrand,
        tname,
        tlen,
        tstart,
        tend,
        matches,
        total_bp,
        qual,
        kwattr,
    ):
        self.qname = qname
        self.qlen = int(qlen)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.samestrand = samestrand
        self.tname = tname
        self.tlen = int(tlen)
        self.tstart = int(tstart)
        self.tend = int(tend)
        self.matches = int(matches)
        self.total_bp = int(total_bp)
        self.qual = int(qual)
        self.kwattr = kwattr
        self.gen_kw_attr()

    def gen_kw_attr(self):
        kwattr_dict = {kw.split(":")[0]: kw.split(":")[-1] for kw in self.kwattr}
        for key in kwattr_dict:
            self.__dict__[key] = kwattr_dict[key]


class Primer:
    def __init__(self, contig, start, end, name, score, strand):
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = int(score)
        self.strand = strand
        self.type = "primary"
        self.amp_no, self.pos = self.get_name_infos()

    def get_name_infos(self):
        ls = self.name.split("_")
        if len(ls) == 4:
            self.type = "alt"
            self.alt_name = ls[3]
        return ls[1], ls[2]


class Amp:
    def __init__(self, amp_no, primers):
        self.name = int(amp_no)
        self.primers = primers
        self.start = min([primer.start for primer in primers])
        self.end = max([primer.end for primer in primers])
        self.max_len = self.end - self.start
        self.fwp_boundary = max(
            [prim for prim in primers if prim.pos == "LEFT"], key=lambda x: x.end
        ).end
        self.revp_boundary = min(
            [prim for prim in primers if prim.pos == "RIGHT"], key=lambda x: x.start
        ).start
        self.read_names = list()
        self.reads = list()

    def random_sample(self, cutoff):
        if len(self.read_names) > cutoff:
            self.selected = random.choices(self.read_names, k=cutoff)
        else:
            self.selected = self.read_names


def create_primer_objs(primer_bed):
    with open(primer_bed, "r") as bed:
        primers = list()
        for line in bed:
            prim = Primer(*line.strip().split("\t"))
            primers.append(prim)
    return sorted(primers, key=lambda x: x.end)


def generate_amps(primers):
    amp_nums = set([primer.amp_no for primer in primers])
    amps = list()
    for num in amp_nums:
        ao = Amp(num, [primer for primer in primers if primer.amp_no == num])
        print(ao.name, ao.max_len)
        amps.append(ao)
    return sorted(amps, key=lambda x: x.name)


def filter_read_mappings(mappings):
    mappings = [m for m in mappings if 300 < m.qlen < 600]
    mappings = [m for m in mappings if m.qual == 60]
    # mappings = [m for m in mappings if m.tp == "P"]
    return mappings


def create_read_mappings(mm2_paf):
    with open(mm2_paf, "r") as paf:
        map_dct = defaultdict(list)
        for line in paf:
            if len(line.strip()) > 0:
                ls = line.strip().split("\t")
                mapping = Mapping(*ls[:12], ls[12:])
                map_dct[mapping.qname].append(mapping)
        mult_maps = {n: ml for n, ml in map_dct.items() if len(ml) > 1}
        mappings = [m for k, l in map_dct.items() for m in l]
        mappings = filter_read_mappings(mappings)
        mono_mappings = [m for m in mappings if m.qname not in mult_maps]
        dual_mappings = {k: v for k, v in mult_maps.items() if len(v) == 2}
        incl_max = [
            max(dual_mappings[mname], key=lambda x: x.matches)
            for mname in dual_mappings
        ]
        incl_max = filter_read_mappings(incl_max)
        mono_mappings.extend(incl_max)
        mappings = mono_mappings
    return sorted(mappings, key=lambda x: x.tend)


def bin_mappings(amp_bins, mappings):
    binned = list()
    na = list()
    while len(amp_bins) > 0:
        if len(mappings) > 0:
            if mappings[0].tend <= amp_bins[0].end + 5:
                if mappings[0].tstart >= amp_bins[0].start - 5:
                    amp_bins[0].read_names.append(mappings[0].qname)
                    mappings.pop(0)
                else:
                    na.append(mappings[0].qname)
                    mappings.pop(0)
            else:
                binned.append(amp_bins[0])
                amp_bins.pop(0)
        else:
            break

    for bin in binned:
        bin.random_sample(200)
        print(bin.name, len(bin.read_names), "selected:", len(bin.selected))
    print("na", len(na))

    return binned


def write_capped_reads(binned, reads, fa_out):
    # print("Writing json")
    # bins_dct = {amp.name:amp.read_names for amp in binned}
    # with open(js_out, "w") as js:
    #     json.dump(bins_dct, js)

    print("Writing fasta")
    all_picks = ["@" + name for amp in binned for name in amp.selected]
    with open(reads, "r") as fq, open(fa_out, "w") as fa:
        for line in fq:
            if line.startswith("@"):
                readname = line.split(" ")[0]
                if readname in all_picks:
                    readname = readname.replace("@", ">")
                    fa.write(readname + "\n")
                    fa.write(next(fq))


if __name__ == "__main__":
    import sys

    # mm2_paf = sys.argv[1]
    # primer_bed = sys.argv[2]
    # reads = sys.argv[3]
    # fa_out = reads + "_capped"
    # js_out = reads + "_capped.json"

    primer_bed = snakemake.input[0]
    mm2_paf = snakemake.input[1]
    reads = snakemake.input[2]
    fa_out = snakemake.output[0]

    primers = create_primer_objs(primer_bed)
    amps = generate_amps(primers)
    mappings = create_read_mappings(mm2_paf)
    binned = bin_mappings(amps, mappings)
    write_capped_reads(binned, reads, fa_out)
