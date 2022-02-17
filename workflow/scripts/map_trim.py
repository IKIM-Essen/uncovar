# import gzip
from collections import defaultdict


class Read:
    def __init__(self, header, seq):
        self.header = header
        self.name = self.header.split(">")[1]
        self.seq = seq
        # self.prim_clipped_seq = ""

    def clip_primers(self, fwp_boundary, revp_boundary, mapping):
        qlen, qstart, qend = mapping.qlen, mapping.qstart, mapping.qend
        tstart, tend = mapping.tstart, mapping.tend

        clip_left = 0
        if tstart <= fwp_boundary:
            add_clip = fwp_boundary - tstart
            clip_left = qstart + add_clip
        else:
            ldiff = tstart - fwp_boundary
            if qstart >= ldiff:
                clip_left = qstart - ldiff

        clip_right = qlen
        if tend >= revp_boundary:
            sub_clip = tend - revp_boundary
            clip_right = qend - sub_clip
        else:
            rdiff = revp_boundary - tend
            if qlen - qend >= rdiff:
                clip_right = qend + rdiff

        self.seq = self.seq[clip_left:clip_right]
        print(clip_left, qlen - clip_right)
        return clip_left, qlen - clip_right


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
        # self.read_names = list()
        self.mappings = dict()
        self.reads = dict()

    def primer_clipping_all(self):
        clip_ct_left = 0
        clip_ct_right = 0
        for read in self.reads:
            try:
                mapping = self.mappings[read]
                left, right = self.reads[read].clip_primers(
                    self.fwp_boundary, self.revp_boundary, mapping
                )
                clip_ct_left += left
                clip_ct_right += right
            except KeyError as e:
                print(f"KeyError in primer_clipping_all: {e}")
        return clip_ct_left, clip_ct_right


def create_primer_objs(primer_bed):
    with open(primer_bed, "r") as bed:
        primers = list()
        for line in bed:
            prim = Primer(*line.strip().split("\t"))
            primers.append(prim)
    return sorted(primers, key=lambda x: x.end)


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


def filter_read_mappings(mappings):
    mappings = [m for m in mappings if 300 < m.qlen < 600]
    mappings = [m for m in mappings if m.qual == 60]
    mappings = [m for m in mappings if m.tp == "P"]
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
    return sorted(mappings, key=lambda x: (x.tend, x.tstart))


def bin_mappings(amp_bins, mappings):
    binned = list()
    na = list()
    # print(amp_bins)
    while len(amp_bins) > 0:
        if len(mappings) > 0:
            if mappings[0].tend <= amp_bins[0].end + 5:
                if mappings[0].tstart >= amp_bins[0].start - 5:
                    amp_bins[0].mappings[mappings[0].qname] = mappings[0]
                    mappings.pop(0)
                else:
                    na.append(mappings[0].qname)
                    mappings.pop(0)
            else:
                binned.append(amp_bins[0])
                amp_bins.pop(0)
        else:
            break

    # for bin in binned:
    #     print(bin.name, "reads:", len(bin.reads))
    # print("na", len(na))

    return binned


def load_reads(read_fasta, amp_bins):

    reads = dict()
    with open(read_fasta, "r") as rfa:
        for line in rfa:
            if line.startswith(">"):
                # print(line)
                header = line.strip().split(" ")[0]
                seq = next(rfa)
                read = Read(header, seq.strip())
                reads[read.name] = read

    # print(reads)

    for amp in amp_bins:
        print("amp.mappings", len(amp.mappings))
        amp.reads = {k: v for k, v in reads.items() if k in amp.mappings}
        print("amp.reads", len(amp.reads))
        # print(amp.reads)

    return amp_bins


def clip_and_write_out(amp_bins, clipped_out):
    with open(clipped_out, "w") as out:
        clip_ct = {"left": 0, "right": 0}
        for amp in amp_bins:
            left, right = amp.primer_clipping_all()
            clip_ct["left"] += left
            clip_ct["right"] += right
            for read in amp.reads:
                out.write(amp.reads[read].header + "\n")
                out.write(amp.reads[read].seq + "\n")
    print(
        f"{clip_ct['left']} bases were clipped from the left/start of reads and "
        f"{clip_ct['right']} bases were clipped from the right/end of reads"
    )


if __name__ == "__main__":
    import sys

    mm2_paf = sys.argv[1]
    primer_bed = sys.argv[2]
    reads = sys.argv[3]
    clipped_out = reads + "_primer_clipped"

    # primer_bed = snakemake.input[0]
    # mm2_paf = snakemake.input[1]
    # reads = snakemake.input[2]

    primers = create_primer_objs(primer_bed)
    amps = generate_amps(primers)
    mappings = create_read_mappings(mm2_paf)
    amps_bin_maps = bin_mappings(amps, mappings)
    amps_bin_reads = load_reads(reads, amps_bin_maps)
    clip_and_write_out(amps_bin_reads, clipped_out)
