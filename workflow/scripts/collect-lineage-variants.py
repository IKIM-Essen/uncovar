# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys
from collections import defaultdict, namedtuple
from enum import Enum
from itertools import product

sys.stderr = open(snakemake.log[0], "w")

import gffutils
import numpy as np
import requests
from dnachisel.biotools import get_backtranslation_table, translate
from pysam import FastaFile, VariantFile, VariantHeader, VariantRecord
from requests.models import ContentDecodingError

covariants_data = requests.get(
    "https://raw.githubusercontent.com/hodcroftlab/covariants/master/web/public/data/clusters.json"
).json()
translate_aa = get_backtranslation_table("Standard")
gff = gffutils.create_db(snakemake.input.annotation, dbfn=":memory:")
gene_start = {gene["gene_name"][0]: gene.start for gene in gff.features_of_type("gene")}
gene_end = {gene["gene_name"][0]: gene.end for gene in gff.features_of_type("gene")}


def aa_to_dna(aa_seq):
    return (
        "".join(combination)
        for combination in product(*[translate_aa[aa] for aa in aa_seq])
    )


def codon_equivalence_class(dna_seq):
    aa = translate(dna_seq)
    return aa_to_dna(aa)


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

    def __repr__(self):
        return repr(self.signature())


class NonSynonymousVariant(SynonymousVariant):
    def __init__(self, left, pos, right, gene):
        super().__init__(left, pos, right)
        self.gene = gene

    def __eq__(self, other):
        return super().__eq__(other) and self.gene == other.gene

    def __hash__(self):
        return hash((self.left, self.pos, self.right, self.gene))

    def is_same_feature(self, other):
        return self.gene == other.gene

    def genome_pos(self):
        return gene_start[self.gene] - 1 + (self.pos - 1) * 3

    def signature(self):
        return f"{self.gene}:{self.left}{self.pos}{self.right}"

    def is_in_first_codon(self):
        return self.pos == 1

    def is_in_last_codon(self):
        aa_len = gene_end[self.gene] - gene_start[self.gene] / 3
        assert self.pos <= aa_len
        return self.pos == aa_len


with FastaFile(snakemake.input.reference) as infasta:
    assert infasta.nreferences == 1
    contig = infasta.references[0]
    ref_len = infasta.lengths[0]
    header = VariantHeader()
    header.add_line(f"##contig=<ID={contig},length={ref_len}")
    header.add_line(
        '##INFO=<ID=SIGNATURES,Number=.,Type=String,Description="Variant signature as obtained from covariants.org">'
    )
    header.add_line(
        '##INFO=<ID=LINEAGES,Number=.,Type=String,Description="Lineages having this variant">'
    )

    known_non_synonymous_variants = defaultdict(set)
    for lineage_entry in covariants_data["clusters"]:
        if (
            "mutations" in lineage_entry
            and "nonsynonymous" in lineage_entry["mutations"]
        ):
            for variant in lineage_entry["mutations"]["nonsynonymous"]:
                variant = NonSynonymousVariant(**variant)
                if variant.gene in gene_start:
                    known_non_synonymous_variants[variant].add(
                        lineage_entry["build_name"]
                    )
                else:
                    print(
                        f"Skipping variant at {variant.gene} because gene is not in given GFF annotation.",
                        file=sys.stderr,
                    )

    known_synonymous_variants = defaultdict(set)
    for lineage_entry in covariants_data["clusters"]:
        if "mutations" in lineage_entry and "synonymous" in lineage_entry["mutations"]:
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
                yield from filtered_variants
            else:

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
                    if len(lineage_matrix) > 0:
                        lineage_matrix = np.unique(lineage_matrix, axis=1)
                    # Step 4: iterate over combinations
                    batch = np.array(batch)
                    batch_lineages = np.array(batch_lineages)
                    for variant_combination in lineage_matrix.T:
                        # select variants and lineages
                        variants = batch[variant_combination]
                        lineages = set.intersection(
                            *batch_lineages[variant_combination]
                        )
                        # yield them in consecutive inner batches
                        last_pos = None
                        inner_batch_start = 0
                        for i, variant in enumerate(variants):
                            if last_pos is not None and variant.pos != last_pos + 1:
                                # yield inner batch
                                yield variants[inner_batch_start:i], lineages
                                inner_batch_start = i
                            last_pos = variant.pos
                        yield variants[inner_batch_start:], lineages

                batch = []
                batch_lineages = []
                for variant, lineages in filtered_variants:
                    if not batch or (
                        variant.pos == batch[-1].pos + 1
                        and variant.is_same_feature(batch[-1])
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
            record.pos = pos + 1  # pysam expects 1-based positions here

            record.info["LINEAGES"] = ",".join(lineages)

            record.info["SIGNATURES"] = ",".join(
                variant.signature() for variant in variants
            )

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
            pos = variant.genome_pos()
            write_record(pos, variant.left, variant.right, lineages, [variant])

        for variants, lineages in get_variants(
            known_non_synonymous_variants, VariantType.Ins
        ):
            pos = variants[0].genome_pos()
            assert not variants[
                0
            ].is_in_first_codon(), "unsupported insertion: is in first codon of protein"
            assert not variants[
                -1
            ].is_in_last_codon(), "unsupported insertion: is in last codon of protein"

            # METHOD: add an unchanged codon before and after the actual variant
            ref_allele = infasta.fetch(reference=contig, start=pos - 3, end=pos + 3)
            pre_codons = codon_equivalence_class(
                infasta.fetch(reference=contig, start=pos - 3, end=pos)
            )
            post_codons = codon_equivalence_class(
                infasta.fetch(reference=contig, start=pos, end=pos + 3)
            )
            for pre_codon, post_codon in product(pre_codons, post_codons):
                for ins_seq in aa_to_dna(
                    "".join(variant.right for variant in variants)
                ):
                    alt_allele = pre_codon + ins_seq + post_codon
                    write_record(pos - 3, ref_allele, alt_allele, lineages, variants)

        for variants, lineages in get_variants(
            known_non_synonymous_variants, VariantType.Del
        ):
            variant = variants[0]
            pos = variants[0].genome_pos()
            del_len = len(variants) * 3

            assert not variants[
                0
            ].is_in_first_codon(), "unsupported deletion: is in first codon of protein"
            assert not variants[
                -1
            ].is_in_last_codon(), "unsupported deletion: is in last codon of protein"

            # METHOD: add an unchanged codon before and after the actual variant
            # in order to capture ambiguity in the alignment
            # before the potential deletion
            pre_codons = codon_equivalence_class(
                infasta.fetch(reference=contig, start=pos - 3, end=pos)
            )
            post_codons = codon_equivalence_class(
                infasta.fetch(
                    reference=contig, start=pos + del_len, end=pos + del_len + 3
                )
            )
            # ref allele including the unchanged codons
            ref_allele = infasta.fetch(
                reference=contig, start=pos - 3, end=pos + del_len + 3
            )
            for pre_codon, post_codon in product(pre_codons, post_codons):
                alt_allele = pre_codon + post_codon
                write_record(pos - 3, ref_allele, alt_allele, lineages, variants)

        for variant, lineages in get_variants(
            known_non_synonymous_variants, VariantType.Subst, merge=False
        ):
            pos = variant.genome_pos()

            ref_allele = infasta.fetch(reference=contig, start=pos, end=pos + 3)
            for alt_allele in aa_to_dna(variant.right):

                write_record(pos, ref_allele, alt_allele, lineages, [variant])
