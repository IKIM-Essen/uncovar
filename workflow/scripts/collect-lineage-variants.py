from collections import namedtuple, defaultdict
from enum import Enum

import numpy as np
import gffutils
import requests
from requests.models import ContentDecodingError
from pysam import VariantFile, VariantHeader, VariantRecord, FastaFile
from dnachisel.biotools import get_backtranslation_table

# TODO: Credits to covariants.org
covariants_data = requests.get(
    "https://github.com/hodcroftlab/covariants/blob/master/web/data/clusters.json"
).json()
translate_aa = get_backtranslation_table("Standard")
gff = gffutils.create_db(snakemake.input.annotation, dbfn=":memory:")


def aa_to_dna(aa_seq):
    return "".join(translate_aa(aa) for aa in aa_seq)


class VariantType(Enum):
    Ins = 1
    Del = 2
    Subst = 3


class SynonymousVariant:
    def __init__(self, left, pos, right):
        self.left = left
        self.pos = pos
        self.right = right

    def __eq__(self, other):
        return (
            self.left == other.left
            and self.right == other.right
            and self.pos == other.pos
        )

    def __hash__(self):
        return hash((self.left, self.pos, self.right))

    def __lt__(self, other):
        return self.pos < other.pos

    def is_same_feature(self, other):
        return True

    def variant_type(self):
        if self.left == "-":
            return VariantType.Ins
        elif self.right == "-":
            return VariantType.Del
        else:
            return VariantType.Subst

    def genome_pos(self):
        return self.pos - 1

    def signature(self):
        return f"{self.left}{self.pos}{self.right}"


class NonSynonymousVariant(SynonymousVariant):
    def __init__(self, left, pos, right, gene):
        super().__init__(self, left, pos, right)
        self.gene = gene

    def __eq__(self, other):
        return super().__eq__(self, other) and self.gene == other.gene

    def __hash__(self):
        return hash((self.left, self.pos, self.right, self.gene))

    def is_same_feature(self, other):
        return self.gene == other.gene

    def genome_pos(self):
        return gff[self.gene].start - 1 + (self.pos - 1) * 3

    def signature(self):
        return f"{self.gene}:{self.left}{self.pos}{self.right}"


with FastaFile(snakemake.input.reference) as infasta:
    assert infasta.nreferences == 1
    contig = infasta.references[0]
    header = VariantHeader()
    header.add_line(f"contig=<ID={contig},length={infasta.lengths[0]}")
    header.add_line(
        'INFO=<ID=SIGNATURES,Number=.,Type=String,Description="Variant signature as obtained from covariants.org">'
    )
    header.add_line(
        'INFO=<ID=LINEAGES,Number=.,Type=String,Description="Lineages having this variant">'
    )

    NonSynonymousVariant = namedtuple(
        "NonSynonymousVariant", ["gene", "left", "pos", "right"]
    )
    SynonymousVariant = namedtuple("SynonymousVariant", ["left", "pos", "right"])

    known_non_synonymous_variants = defaultdict(set)
    for lineage_entry in covariants_data["clusters"]:
        for variant in lineage_entry["mutations"]["nonsynonymous"]:
            known_non_synonymous_variants[NonSynonymousVariant(**variant)].add(
                lineage_entry["build_name"]
            )

    known_synonymous_variants = defaultdict(set)
    for lineage_entry in covariants_data["clusters"]:
        for variant in lineage_entry["mutations"]["synonymous"]:
            known_synonymous_variants[SynonymousVariant(**variant)].add(
                lineage_entry["build_name"]
            )

    with VariantFile(snakemake.output[0], "wb", header=header) as outvcf:

        def get_variants(all_variants, variant_type, merge=True):
            filtered_variants = sorted(
                filter(
                    lambda item: item[0].variant_type() == variant_type,
                    all_variants.items(),
                )
            )

            if not merge:
                return filtered_variants

            def process_batch(batch, batch_lineages):
                # Step 1: collect all visited lineages in batch
                all_lineages = np.array(
                    list(
                        set(
                            lineage
                            for lineages in batch_lineages
                            for lineage in lineages
                        )
                    )
                )
                # Step 2: build matrix of variants vs lineages (columns mark combinations of variants that can be merged)
                lineage_matrix = np.array(
                    [
                        [(lineage in lineages) for lineage in all_lineages]
                        for lineages in batch_lineages
                    ]
                )
                # Step 3: remove duplicate columns
                lineage_matrix = np.unique(lineage_matrix, axis=1)
                # Step 4: iterate over combinations
                batch = np.array(batch)
                batch_lineages = np.array(batch_lineages)
                for variant_combination in lineage_matrix.T:
                    # select variants and lineages
                    variants = batch[variant_combination]
                    lineages = set.intersection(*batch_lineages[variant_combination])
                    # yield them in consecutive inner batches
                    last_pos = None
                    inner_batch_start = 0
                    for i, variant in enumerate(variants):
                        if variant.pos != last_pos + 1:
                            # yield inner batch
                            yield variants[inner_batch_start:i], lineages
                            inner_batch_start = i
                        last_pos = variant.pos
                    yield variants[inner_batch_start:], lineages

            batch = []
            batch_lineages = []
            for variant, lineages in filtered_variants:
                if variant.pos == batch[-1].pos + 1 and variant.is_same_feature(
                    batch[-1]
                ):
                    batch.append(variant)
                    batch_lineages.append(lineages)
                else:
                    # yield and remove the last batch
                    yield from process_batch(batch, batch_lineages)
                    # clear and start with new batch
                    batch = [variant]
                    batch_lineages = [lineages]

            yield from process_batch(batch, batch_lineages)

        def write_record(pos, ref_allele, alt_allele, lineages, variants):
            record = outvcf.new_record()
            record.contig = contig
            record.alleles = (ref_allele, alt_allele)
            record.pos = pos

            record.info["LINEAGES"] = lineages

            record.info["SIGNATURES"] = [variant.signature() for variant in variants]

            outvcf.write(record)

        for variants, lineages in get_variants(
            known_synonymous_variants, VariantType.Ins
        ):
            pos = variants[0].genome_pos() - 1
            ref_allele = infasta.fetch(reference=contig, start=pos, end=pos + 1)
            alt_allele = ref_allele + "".join(variant.right for variant in variants)
            write_record(pos, ref_allele, alt_allele, lineages, variants)

        for variants, lineages in get_variants(
            known_synonymous_variants, VariantType.Del
        ):
            pos = variants[0].genome_pos() - 1
            alt_allele = infasta.fetch(reference=contig, start=pos, end=pos + 1)
            ref_allele = alt_allele + "".join(variant.left for variant in variants)
            write_record(pos, ref_allele, alt_allele, lineages, variants)

        for variant, lineages in get_variants(
            known_synonymous_variants, VariantType.Subst, merge=False
        ):
            pos = variants[0].genome_pos()
            write_record(pos, variant.left, variant.right, lineages, [variant])

        for variants, lineages in get_variants(
            known_non_synonymous_variants, VariantType.Ins
        ):
            pos = variants[0].genome_pos() - 1
            ref_allele = infasta.fetch(reference=contig, start=pos, end=pos + 1)
            for ins_seq in aa_to_dna("".join(variant.right for variant in variants)):
                alt_allele = ref_allele + ins_seq
                write_record(pos, ref_allele, alt_allele, lineages, variants)

        for variants, lineages in get_variants(
            known_non_synonymous_variants, VariantType.Del
        ):
            pos = variants[0].genome_pos() - 1
            alt_allele = infasta.fetch(reference=contig, start=pos, end=pos + 1)
            ref_allele = alt_allele + aa_to_dna(
                "".join(variant.left for variant in variants)
            )
            write_record(pos, ref_allele, alt_allele, lineages, variants)

        for variant, lineages in get_variants(
            known_non_synonymous_variants, VariantType.Subst, merge=False
        ):
            pos = variants.genome_pos()

            ref_allele = infasta.fetch(reference=contig, start=pos, end=pos + 3)
            for alt_allele in aa_to_dna(variant.right):
                write_record(pos, ref_allele, alt_allele, lineages, [variant])
