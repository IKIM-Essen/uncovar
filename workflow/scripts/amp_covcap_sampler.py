import random


class Mapping:
    def __init__(self, qname, qlen, qstart, qend, samestrand, tname, tlen, tstart, tend, matches, total_bp, qual, kwattr):
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
            self.key = kwattr_dict[key]


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
        self.reads = list()
        

    def random_sample(self, cutoff):
        if len(self.reads) > cutoff:
            self.selected = random.choices(self.reads, k=cutoff)
        else:
            self.selected = self.reads
            


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
        amps.append(ao)
    return sorted(amps, key=lambda x: x.name)

        
def create_read_mappings(mm2_paf):
    with open(mm2_paf, "r") as paf:
        mappings = list()
        for line in paf:
            if len(line.strip()) > 0:
                ls = line.strip().split("\t")
                mapping = Mapping(*ls[:12], ls[12:])
                mappings.append(mapping)
    return sorted(mappings, key=lambda x: x.tend)


def bin_mappings(amp_bins, mappings):
    binned = list()
    na = list()
    while len(mappings) > 0:
        if mappings[0].tend <= amp_bins[0].end + 5:
            if mappings[0].tstart >= amp_bins[0].start - 5:
                amp_bins[0].reads.append(mappings[0].qname)
                mappings.pop(0)
            else:
                na.append(mappings[0].qname)
                mappings.pop(0)
        else:
            binned.append(amp_bins[0])
            amp_bins.pop(0)

    for bin in binned:
        bin.random_sample(10)
        print(bin.name, len(bin.reads), "selected:", len(bin.selected))
    print("na", len(na))

    return binned

def write_capped_reads(binned, reads, out):
    all_picks = ["@"+ name for amp in binned for name in amp.selected]
    print(all_picks)
    with open(reads, "r") as fq, open(out, "w") as fa:
        for line in fq:
            if line.startswith("@"):
                readname = line.split(" ")[0]
                print(readname)
                if readname in all_picks:
                    readname = readname.replace("@", ">")
                    fa.write(readname + "\n")
                    fa.write(next(fq))


if __name__ == "__main__":
    import sys

    mm2_paf = sys.argv[1]
    primer_bed = sys.argv[2]
    reads = sys.argv[3]
    out = sys.argv[3] + "_capped"

    # primer_bed = snakemake.output[0]
    # mm2_paf = snakemake.input[1]
    # reads = snakemake.input[2]


    primers = create_primer_objs(primer_bed)
    amps = generate_amps(primers)
    mappings = create_read_mappings(mm2_paf)
    binned = bin_mappings(amps, mappings)
    write_capped_reads(binned, reads, out)
        