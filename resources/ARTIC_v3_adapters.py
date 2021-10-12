"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains the class and sequences for known adapters used in Oxford Nanopore library
preparation kits.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""


class Adapter(object):

    def __init__(self, name, start_sequence=None, end_sequence=None, both_ends_sequence=None):
        self.name = name
        self.start_sequence = start_sequence if start_sequence else []
        self.end_sequence = end_sequence if end_sequence else []
        if both_ends_sequence:
            self.start_sequence = both_ends_sequence
            self.end_sequence = both_ends_sequence
        self.best_start_score, self.best_end_score = 0.0, 0.0

    def best_start_or_end_score(self):
        return max(self.best_start_score, self.best_end_score)

    def is_barcode(self):
        return self.name.startswith('Barcode ')

    def barcode_direction(self):
        if '_rev' in self.start_sequence[0]:
            return 'reverse'
        else:
            return 'forward'

    def get_barcode_name(self):
        """
        Gets the barcode name for the output files. We want a concise name, so it looks at all
        options and chooses the shortest.
        """
        possible_names = [self.name]
        if self.start_sequence:
            possible_names.append(self.start_sequence[0])
        if self.end_sequence:
            possible_names.append(self.end_sequence[0])
        barcode_name = sorted(possible_names, key=lambda x: len(x))[0]
        return barcode_name.replace(' ', '_')


# INSTRUCTIONS FOR ADDING CUSTOM ADAPTERS
# ---------------------------------------
# If you need Porechop to remove adapters that aren't included, you can add your own my modifying
# the ADAPTERS list below.
#
# Here is the format for a normal adapter:
#     Adapter('Adapter_set_name',
#             start_sequence=('Start_adapter_name', 'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'),
#             end_sequence=('End_adapter_name', 'AACCGGTTAACCGGTTAACCGGTTAACCGGTT'))
#
# You can exclude start_sequence and end_sequence as appropriate.
#
# If you have custom Barcodes, make sure that the adapter set name starts with 'Barcode '. Also,
# remove the existing barcode sequences from this file to avoid conflicts:
#     Adapter('Barcode 1',
#             start_sequence=('Barcode_1_start', 'AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT'),
#             end_sequence=('Barcode_1_end', 'AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT')),
#     Adapter('Barcode 2',
#             start_sequence=('Barcode_2_start', 'TTTTTTTTGGGGGGGGCCCCCCCCAAAAAAAA'),
#             end_sequence=('Barcode_2_end', 'TTTTTTTTGGGGGGGGCCCCCCCCAAAAAAAA'))


ADAPTERS = [
            Adapter('nCoV-2019_1_LEFT',
                start_sequence=('nCoV-2019_1_LEFT_fw', 'ACCAACCAACTTTCGATCTCTTGT'),
                end_sequence=('nCoV-2019_1_LEFT_revcomp', 'ACAAGAGATCGAAAGTTGGTTGGT')),

            Adapter('nCoV-2019_1_RIGHT',
                start_sequence=('nCoV-2019_1_RIGHT_fw', 'CATCTTTAAGATGTTGACGTGCCTC'),
                end_sequence=('nCoV-2019_1_RIGHT_revcomp', 'GAGGCACGTCAACATCTTAAAGATG')),

            Adapter('nCoV-2019_2_LEFT',
                start_sequence=('nCoV-2019_2_LEFT_fw', 'CTGTTTTACAGGTTCGCGACGT'),
                end_sequence=('nCoV-2019_2_LEFT_revcomp', 'ACGTCGCGAACCTGTAAAACAG')),

            Adapter('nCoV-2019_2_RIGHT',
                start_sequence=('nCoV-2019_2_RIGHT_fw', 'TAAGGATCAGTGCCAAGCTCGT'),
                end_sequence=('nCoV-2019_2_RIGHT_revcomp', 'ACGAGCTTGGCACTGATCCTTA')),

            Adapter('nCoV-2019_3_LEFT',
                start_sequence=('nCoV-2019_3_LEFT_fw', 'CGGTAATAAAGGAGCTGGTGGC'),
                end_sequence=('nCoV-2019_3_LEFT_revcomp', 'GCCACCAGCTCCTTTATTACCG')),

            Adapter('nCoV-2019_3_RIGHT',
                start_sequence=('nCoV-2019_3_RIGHT_fw', 'AAGGTGTCTGCAATTCATAGCTCT'),
                end_sequence=('nCoV-2019_3_RIGHT_revcomp', 'AGAGCTATGAATTGCAGACACCTT')),

            Adapter('nCoV-2019_4_LEFT',
                start_sequence=('nCoV-2019_4_LEFT_fw', 'GGTGTATACTGCTGCCGTGAAC'),
                end_sequence=('nCoV-2019_4_LEFT_revcomp', 'GTTCACGGCAGCAGTATACACC')),

            Adapter('nCoV-2019_4_RIGHT',
                start_sequence=('nCoV-2019_4_RIGHT_fw', 'CACAAGTAGTGGCACCTTCTTTAGT'),
                end_sequence=('nCoV-2019_4_RIGHT_revcomp', 'ACTAAAGAAGGTGCCACTACTTGTG')),

            Adapter('nCoV-2019_5_LEFT',
                start_sequence=('nCoV-2019_5_LEFT_fw', 'TGGTGAAACTTCATGGCAGACG'),
                end_sequence=('nCoV-2019_5_LEFT_revcomp', 'CGTCTGCCATGAAGTTTCACCA')),

            Adapter('nCoV-2019_5_RIGHT',
                start_sequence=('nCoV-2019_5_RIGHT_fw', 'ATTGATGTTGACTTTCTCTTTTTGGAGT'),
                end_sequence=('nCoV-2019_5_RIGHT_revcomp', 'ACTCCAAAAAGAGAAAGTCAACATCAAT')),

            Adapter('nCoV-2019_6_LEFT',
                start_sequence=('nCoV-2019_6_LEFT_fw', 'GGTGTTGTTGGAGAAGGTTCCG'),
                end_sequence=('nCoV-2019_6_LEFT_revcomp', 'CGGAACCTTCTCCAACAACACC')),

            Adapter('nCoV-2019_6_RIGHT',
                start_sequence=('nCoV-2019_6_RIGHT_fw', 'TAGCGGCCTTCTGTAAAACACG'),
                end_sequence=('nCoV-2019_6_RIGHT_revcomp', 'CGTGTTTTACAGAAGGCCGCTA')),

            Adapter('nCoV-2019_7_LEFT',
                start_sequence=('nCoV-2019_7_LEFT_fw', 'ATCAGAGGCTGCTCGTGTTGTA'),
                end_sequence=('nCoV-2019_7_LEFT_revcomp', 'TACAACACGAGCAGCCTCTGAT')),

            Adapter('nCoV-2019_7_LEFT_alt0',
                start_sequence=('nCoV-2019_7_LEFT_alt0_fw', 'CATTTGCATCAGAGGCTGCTCG'),
                end_sequence=('nCoV-2019_7_LEFT_alt0_revcomp', 'CGAGCAGCCTCTGATGCAAATG')),

            Adapter('nCoV-2019_7_RIGHT',
                start_sequence=('nCoV-2019_7_RIGHT_fw', 'TGCACAGGTGACAATTTGTCCA'),
                end_sequence=('nCoV-2019_7_RIGHT_revcomp', 'TGGACAAATTGTCACCTGTGCA')),

            Adapter('nCoV-2019_7_RIGHT_alt5',
                start_sequence=('nCoV-2019_7_RIGHT_alt5_fw', 'AGGTGACAATTTGTCCACCGAC'),
                end_sequence=('nCoV-2019_7_RIGHT_alt5_revcomp', 'GTCGGTGGACAAATTGTCACCT')),

            Adapter('nCoV-2019_8_LEFT',
                start_sequence=('nCoV-2019_8_LEFT_fw', 'AGAGTTTCTTAGAGACGGTTGGGA'),
                end_sequence=('nCoV-2019_8_LEFT_revcomp', 'TCCCAACCGTCTCTAAGAAACTCT')),

            Adapter('nCoV-2019_8_RIGHT',
                start_sequence=('nCoV-2019_8_RIGHT_fw', 'GCTTCAACAGCTTCACTAGTAGGT'),
                end_sequence=('nCoV-2019_8_RIGHT_revcomp', 'ACCTACTAGTGAAGCTGTTGAAGC')),

            Adapter('nCoV-2019_9_LEFT',
                start_sequence=('nCoV-2019_9_LEFT_fw', 'TCCCACAGAAGTGTTAACAGAGGA'),
                end_sequence=('nCoV-2019_9_LEFT_revcomp', 'TCCTCTGTTAACACTTCTGTGGGA')),

            Adapter('nCoV-2019_9_LEFT_alt4',
                start_sequence=('nCoV-2019_9_LEFT_alt4_fw', 'TTCCCACAGAAGTGTTAACAGAGG'),
                end_sequence=('nCoV-2019_9_LEFT_alt4_revcomp', 'CCTCTGTTAACACTTCTGTGGGAA')),

            Adapter('nCoV-2019_9_RIGHT',
                start_sequence=('nCoV-2019_9_RIGHT_fw', 'ATGACAGCATCTGCCACAACAC'),
                end_sequence=('nCoV-2019_9_RIGHT_revcomp', 'GTGTTGTGGCAGATGCTGTCAT')),

            Adapter('nCoV-2019_9_RIGHT_alt2',
                start_sequence=('nCoV-2019_9_RIGHT_alt2_fw', 'GACAGCATCTGCCACAACACAG'),
                end_sequence=('nCoV-2019_9_RIGHT_alt2_revcomp', 'CTGTGTTGTGGCAGATGCTGTC')),

            Adapter('nCoV-2019_10_LEFT',
                start_sequence=('nCoV-2019_10_LEFT_fw', 'TGAGAAGTGCTCTGCCTATACAGT'),
                end_sequence=('nCoV-2019_10_LEFT_revcomp', 'ACTGTATAGGCAGAGCACTTCTCA')),

            Adapter('nCoV-2019_10_RIGHT',
                start_sequence=('nCoV-2019_10_RIGHT_fw', 'TCATCTAACCAATCTTCTTCTTGCTCT'),
                end_sequence=('nCoV-2019_10_RIGHT_revcomp', 'AGAGCAAGAAGAAGATTGGTTAGATGA')),

            Adapter('nCoV-2019_11_LEFT',
                start_sequence=('nCoV-2019_11_LEFT_fw', 'GGAATTTGGTGCCACTTCTGCT'),
                end_sequence=('nCoV-2019_11_LEFT_revcomp', 'AGCAGAAGTGGCACCAAATTCC')),

            Adapter('nCoV-2019_11_RIGHT',
                start_sequence=('nCoV-2019_11_RIGHT_fw', 'TCATCAGATTCAACTTGCATGGCA'),
                end_sequence=('nCoV-2019_11_RIGHT_revcomp', 'TGCCATGCAAGTTGAATCTGATGA')),

            Adapter('nCoV-2019_12_LEFT',
                start_sequence=('nCoV-2019_12_LEFT_fw', 'AAACATGGAGGAGGTGTTGCAG'),
                end_sequence=('nCoV-2019_12_LEFT_revcomp', 'CTGCAACACCTCCTCCATGTTT')),

            Adapter('nCoV-2019_12_RIGHT',
                start_sequence=('nCoV-2019_12_RIGHT_fw', 'TTCACTCTTCATTTCCAAAAAGCTTGA'),
                end_sequence=('nCoV-2019_12_RIGHT_revcomp', 'TCAAGCTTTTTGGAAATGAAGAGTGAA')),

            Adapter('nCoV-2019_13_LEFT',
                start_sequence=('nCoV-2019_13_LEFT_fw', 'TCGCACAAATGTCTACTTAGCTGT'),
                end_sequence=('nCoV-2019_13_LEFT_revcomp', 'ACAGCTAAGTAGACATTTGTGCGA')),

            Adapter('nCoV-2019_13_RIGHT',
                start_sequence=('nCoV-2019_13_RIGHT_fw', 'ACCACAGCAGTTAAAACACCCT'),
                end_sequence=('nCoV-2019_13_RIGHT_revcomp', 'AGGGTGTTTTAACTGCTGTGGT')),

            Adapter('nCoV-2019_14_LEFT',
                start_sequence=('nCoV-2019_14_LEFT_fw', 'CATCCAGATTCTGCCACTCTTGT'),
                end_sequence=('nCoV-2019_14_LEFT_revcomp', 'ACAAGAGTGGCAGAATCTGGATG')),

            Adapter('nCoV-2019_14_LEFT_alt4',
                start_sequence=('nCoV-2019_14_LEFT_alt4_fw', 'TGGCAATCTTCATCCAGATTCTGC'),
                end_sequence=('nCoV-2019_14_LEFT_alt4_revcomp', 'GCAGAATCTGGATGAAGATTGCCA')),

            Adapter('nCoV-2019_14_RIGHT',
                start_sequence=('nCoV-2019_14_RIGHT_fw', 'AGTTTCCACACAGACAGGCATT'),
                end_sequence=('nCoV-2019_14_RIGHT_revcomp', 'AATGCCTGTCTGTGTGGAAACT')),

            Adapter('nCoV-2019_14_RIGHT_alt2',
                start_sequence=('nCoV-2019_14_RIGHT_alt2_fw', 'TGCGTGTTTCTTCTGCATGTGC'),
                end_sequence=('nCoV-2019_14_RIGHT_alt2_revcomp', 'GCACATGCAGAAGAAACACGCA')),

            Adapter('nCoV-2019_15_LEFT',
                start_sequence=('nCoV-2019_15_LEFT_fw', 'ACAGTGCTTAAAAAGTGTAAAAGTGCC'),
                end_sequence=('nCoV-2019_15_LEFT_revcomp', 'GGCACTTTTACACTTTTTAAGCACTGT')),

            Adapter('nCoV-2019_15_LEFT_alt1',
                start_sequence=('nCoV-2019_15_LEFT_alt1_fw', 'AGTGCTTAAAAAGTGTAAAAGTGCCT'),
                end_sequence=('nCoV-2019_15_LEFT_alt1_revcomp', 'AGGCACTTTTACACTTTTTAAGCACT')),

            Adapter('nCoV-2019_15_RIGHT',
                start_sequence=('nCoV-2019_15_RIGHT_fw', 'AACAGAAACTGTAGCTGGCACT'),
                end_sequence=('nCoV-2019_15_RIGHT_revcomp', 'AGTGCCAGCTACAGTTTCTGTT')),

            Adapter('nCoV-2019_15_RIGHT_alt3',
                start_sequence=('nCoV-2019_15_RIGHT_alt3_fw', 'ACTGTAGCTGGCACTTTGAGAGA'),
                end_sequence=('nCoV-2019_15_RIGHT_alt3_revcomp', 'TCTCTCAAAGTGCCAGCTACAGT')),

            Adapter('nCoV-2019_16_LEFT',
                start_sequence=('nCoV-2019_16_LEFT_fw', 'AATTTGGAAGAAGCTGCTCGGT'),
                end_sequence=('nCoV-2019_16_LEFT_revcomp', 'ACCGAGCAGCTTCTTCCAAATT')),

            Adapter('nCoV-2019_16_RIGHT',
                start_sequence=('nCoV-2019_16_RIGHT_fw', 'CACAACTTGCGTGTGGAGGTTA'),
                end_sequence=('nCoV-2019_16_RIGHT_revcomp', 'TAACCTCCACACGCAAGTTGTG')),

            Adapter('nCoV-2019_17_LEFT',
                start_sequence=('nCoV-2019_17_LEFT_fw', 'CTTCTTTCTTTGAGAGAAGTGAGGACT'),
                end_sequence=('nCoV-2019_17_LEFT_revcomp', 'AGTCCTCACTTCTCTCAAAGAAAGAAG')),

            Adapter('nCoV-2019_17_RIGHT',
                start_sequence=('nCoV-2019_17_RIGHT_fw', 'TTTGTTGGAGTGTTAACAATGCAGT'),
                end_sequence=('nCoV-2019_17_RIGHT_revcomp', 'ACTGCATTGTTAACACTCCAACAAA')),

            Adapter('nCoV-2019_18_LEFT',
                start_sequence=('nCoV-2019_18_LEFT_fw', 'TGGAAATACCCACAAGTTAATGGTTTAAC'),
                end_sequence=('nCoV-2019_18_LEFT_revcomp', 'GTTAAACCATTAACTTGTGGGTATTTCCA')),

            Adapter('nCoV-2019_18_LEFT_alt2',
                start_sequence=('nCoV-2019_18_LEFT_alt2_fw', 'ACTTCTATTAAATGGGCAGATAACAACTGT'),
                end_sequence=('nCoV-2019_18_LEFT_alt2_revcomp', 'ACAGTTGTTATCTGCCCATTTAATAGAAGT')),

            Adapter('nCoV-2019_18_RIGHT',
                start_sequence=('nCoV-2019_18_RIGHT_fw', 'AGCTTGTTTACCACACGTACAAGG'),
                end_sequence=('nCoV-2019_18_RIGHT_revcomp', 'CCTTGTACGTGTGGTAAACAAGCT')),

            Adapter('nCoV-2019_18_RIGHT_alt1',
                start_sequence=('nCoV-2019_18_RIGHT_alt1_fw', 'GCTTGTTTACCACACGTACAAGG'),
                end_sequence=('nCoV-2019_18_RIGHT_alt1_revcomp', 'CCTTGTACGTGTGGTAAACAAGC')),

            Adapter('nCoV-2019_19_LEFT',
                start_sequence=('nCoV-2019_19_LEFT_fw', 'GCTGTTATGTACATGGGCACACT'),
                end_sequence=('nCoV-2019_19_LEFT_revcomp', 'AGTGTGCCCATGTACATAACAGC')),

            Adapter('nCoV-2019_19_RIGHT',
                start_sequence=('nCoV-2019_19_RIGHT_fw', 'TGTCCAACTTAGGGTCAATTTCTGT'),
                end_sequence=('nCoV-2019_19_RIGHT_revcomp', 'ACAGAAATTGACCCTAAGTTGGACA')),

            Adapter('nCoV-2019_20_LEFT',
                start_sequence=('nCoV-2019_20_LEFT_fw', 'ACAAAGAAAACAGTTACACAACAACCA'),
                end_sequence=('nCoV-2019_20_LEFT_revcomp', 'TGGTTGTTGTGTAACTGTTTTCTTTGT')),

            Adapter('nCoV-2019_20_RIGHT',
                start_sequence=('nCoV-2019_20_RIGHT_fw', 'ACGTGGCTTTATTAGTTGCATTGTT'),
                end_sequence=('nCoV-2019_20_RIGHT_revcomp', 'AACAATGCAACTAATAAAGCCACGT')),

            Adapter('nCoV-2019_21_LEFT',
                start_sequence=('nCoV-2019_21_LEFT_fw', 'TGGCTATTGATTATAAACACTACACACCC'),
                end_sequence=('nCoV-2019_21_LEFT_revcomp', 'GGGTGTGTAGTGTTTATAATCAATAGCCA')),

            Adapter('nCoV-2019_21_LEFT_alt2',
                start_sequence=('nCoV-2019_21_LEFT_alt2_fw', 'GGCTATTGATTATAAACACTACACACCCT'),
                end_sequence=('nCoV-2019_21_LEFT_alt2_revcomp', 'AGGGTGTGTAGTGTTTATAATCAATAGCC')),

            Adapter('nCoV-2019_21_RIGHT',
                start_sequence=('nCoV-2019_21_RIGHT_fw', 'TAGATCTGTGTGGCCAACCTCT'),
                end_sequence=('nCoV-2019_21_RIGHT_revcomp', 'AGAGGTTGGCCACACAGATCTA')),

            Adapter('nCoV-2019_21_RIGHT_alt0',
                start_sequence=('nCoV-2019_21_RIGHT_alt0_fw', 'GATCTGTGTGGCCAACCTCTTC'),
                end_sequence=('nCoV-2019_21_RIGHT_alt0_revcomp', 'GAAGAGGTTGGCCACACAGATC')),

            Adapter('nCoV-2019_22_LEFT',
                start_sequence=('nCoV-2019_22_LEFT_fw', 'ACTACCGAAGTTGTAGGAGACATTATACT'),
                end_sequence=('nCoV-2019_22_LEFT_revcomp', 'AGTATAATGTCTCCTACAACTTCGGTAGT')),

            Adapter('nCoV-2019_22_RIGHT',
                start_sequence=('nCoV-2019_22_RIGHT_fw', 'ACAGTATTCTTTGCTATAGTAGTCGGC'),
                end_sequence=('nCoV-2019_22_RIGHT_revcomp', 'GCCGACTACTATAGCAAAGAATACTGT')),

            Adapter('nCoV-2019_23_LEFT',
                start_sequence=('nCoV-2019_23_LEFT_fw', 'ACAACTACTAACATAGTTACACGGTGT'),
                end_sequence=('nCoV-2019_23_LEFT_revcomp', 'ACACCGTGTAACTATGTTAGTAGTTGT')),

            Adapter('nCoV-2019_23_RIGHT',
                start_sequence=('nCoV-2019_23_RIGHT_fw', 'ACCAGTACAGTAGGTTGCAATAGTG'),
                end_sequence=('nCoV-2019_23_RIGHT_revcomp', 'CACTATTGCAACCTACTGTACTGGT')),

            Adapter('nCoV-2019_24_LEFT',
                start_sequence=('nCoV-2019_24_LEFT_fw', 'AGGCATGCCTTCTTACTGTACTG'),
                end_sequence=('nCoV-2019_24_LEFT_revcomp', 'CAGTACAGTAAGAAGGCATGCCT')),

            Adapter('nCoV-2019_24_RIGHT',
                start_sequence=('nCoV-2019_24_RIGHT_fw', 'ACATTCTAACCATAGCTGAAATCGGG'),
                end_sequence=('nCoV-2019_24_RIGHT_revcomp', 'CCCGATTTCAGCTATGGTTAGAATGT')),

            Adapter('nCoV-2019_25_LEFT',
                start_sequence=('nCoV-2019_25_LEFT_fw', 'GCAATTGTTTTTCAGCTATTTTGCAGT'),
                end_sequence=('nCoV-2019_25_LEFT_revcomp', 'ACTGCAAAATAGCTGAAAAACAATTGC')),

            Adapter('nCoV-2019_25_RIGHT',
                start_sequence=('nCoV-2019_25_RIGHT_fw', 'ACTGTAGTGACAAGTCTCTCGCA'),
                end_sequence=('nCoV-2019_25_RIGHT_revcomp', 'TGCGAGAGACTTGTCACTACAGT')),

            Adapter('nCoV-2019_26_LEFT',
                start_sequence=('nCoV-2019_26_LEFT_fw', 'TTGTGATACATTCTGTGCTGGTAGT'),
                end_sequence=('nCoV-2019_26_LEFT_revcomp', 'ACTACCAGCACAGAATGTATCACAA')),

            Adapter('nCoV-2019_26_RIGHT',
                start_sequence=('nCoV-2019_26_RIGHT_fw', 'TCCGCACTATCACCAACATCAG'),
                end_sequence=('nCoV-2019_26_RIGHT_revcomp', 'CTGATGTTGGTGATAGTGCGGA')),

            Adapter('nCoV-2019_27_LEFT',
                start_sequence=('nCoV-2019_27_LEFT_fw', 'ACTACAGTCAGCTTATGTGTCAACC'),
                end_sequence=('nCoV-2019_27_LEFT_revcomp', 'GGTTGACACATAAGCTGACTGTAGT')),

            Adapter('nCoV-2019_27_RIGHT',
                start_sequence=('nCoV-2019_27_RIGHT_fw', 'AATACAAGCACCAAGGTCACGG'),
                end_sequence=('nCoV-2019_27_RIGHT_revcomp', 'CCGTGACCTTGGTGCTTGTATT')),

            Adapter('nCoV-2019_28_LEFT',
                start_sequence=('nCoV-2019_28_LEFT_fw', 'ACATAGAAGTTACTGGCGATAGTTGT'),
                end_sequence=('nCoV-2019_28_LEFT_revcomp', 'ACAACTATCGCCAGTAACTTCTATGT')),

            Adapter('nCoV-2019_28_RIGHT',
                start_sequence=('nCoV-2019_28_RIGHT_fw', 'TGTTTAGACATGACATGAACAGGTGT'),
                end_sequence=('nCoV-2019_28_RIGHT_revcomp', 'ACACCTGTTCATGTCATGTCTAAACA')),

            Adapter('nCoV-2019_29_LEFT',
                start_sequence=('nCoV-2019_29_LEFT_fw', 'ACTTGTGTTCCTTTTTGTTGCTGC'),
                end_sequence=('nCoV-2019_29_LEFT_revcomp', 'GCAGCAACAAAAAGGAACACAAGT')),

            Adapter('nCoV-2019_29_RIGHT',
                start_sequence=('nCoV-2019_29_RIGHT_fw', 'AGTGTACTCTATAAGTTTTGATGGTGTGT'),
                end_sequence=('nCoV-2019_29_RIGHT_revcomp', 'ACACACCATCAAAACTTATAGAGTACACT')),

            Adapter('nCoV-2019_30_LEFT',
                start_sequence=('nCoV-2019_30_LEFT_fw', 'GCACAACTAATGGTGACTTTTTGCA'),
                end_sequence=('nCoV-2019_30_LEFT_revcomp', 'TGCAAAAAGTCACCATTAGTTGTGC')),

            Adapter('nCoV-2019_30_RIGHT',
                start_sequence=('nCoV-2019_30_RIGHT_fw', 'ACCACTAGTAGATACACAAACACCAG'),
                end_sequence=('nCoV-2019_30_RIGHT_revcomp', 'CTGGTGTTTGTGTATCTACTAGTGGT')),

            Adapter('nCoV-2019_31_LEFT',
                start_sequence=('nCoV-2019_31_LEFT_fw', 'TTCTGAGTACTGTAGGCACGGC'),
                end_sequence=('nCoV-2019_31_LEFT_revcomp', 'GCCGTGCCTACAGTACTCAGAA')),

            Adapter('nCoV-2019_31_RIGHT',
                start_sequence=('nCoV-2019_31_RIGHT_fw', 'ACAGAATAAACACCAGGTAAGAATGAGT'),
                end_sequence=('nCoV-2019_31_RIGHT_revcomp', 'ACTCATTCTTACCTGGTGTTTATTCTGT')),

            Adapter('nCoV-2019_32_LEFT',
                start_sequence=('nCoV-2019_32_LEFT_fw', 'TGGTGAATACAGTCATGTAGTTGCC'),
                end_sequence=('nCoV-2019_32_LEFT_revcomp', 'GGCAACTACATGACTGTATTCACCA')),

            Adapter('nCoV-2019_32_RIGHT',
                start_sequence=('nCoV-2019_32_RIGHT_fw', 'AGCACATCACTACGCAACTTTAGA'),
                end_sequence=('nCoV-2019_32_RIGHT_revcomp', 'TCTAAAGTTGCGTAGTGATGTGCT')),

            Adapter('nCoV-2019_33_LEFT',
                start_sequence=('nCoV-2019_33_LEFT_fw', 'ACTTTTGAAGAAGCTGCGCTGT'),
                end_sequence=('nCoV-2019_33_LEFT_revcomp', 'ACAGCGCAGCTTCTTCAAAAGT')),

            Adapter('nCoV-2019_33_RIGHT',
                start_sequence=('nCoV-2019_33_RIGHT_fw', 'TGGACAGTAAACTACGTCATCAAGC'),
                end_sequence=('nCoV-2019_33_RIGHT_revcomp', 'GCTTGATGACGTAGTTTACTGTCCA')),

            Adapter('nCoV-2019_34_LEFT',
                start_sequence=('nCoV-2019_34_LEFT_fw', 'TCCCATCTGGTAAAGTTGAGGGT'),
                end_sequence=('nCoV-2019_34_LEFT_revcomp', 'ACCCTCAACTTTACCAGATGGGA')),

            Adapter('nCoV-2019_34_RIGHT',
                start_sequence=('nCoV-2019_34_RIGHT_fw', 'AGTGAAATTGGGCCTCATAGCA'),
                end_sequence=('nCoV-2019_34_RIGHT_revcomp', 'TGCTATGAGGCCCAATTTCACT')),

            Adapter('nCoV-2019_35_LEFT',
                start_sequence=('nCoV-2019_35_LEFT_fw', 'TGTTCGCATTCAACCAGGACAG'),
                end_sequence=('nCoV-2019_35_LEFT_revcomp', 'CTGTCCTGGTTGAATGCGAACA')),

            Adapter('nCoV-2019_35_RIGHT',
                start_sequence=('nCoV-2019_35_RIGHT_fw', 'ACTTCATAGCCACAAGGTTAAAGTCA'),
                end_sequence=('nCoV-2019_35_RIGHT_revcomp', 'TGACTTTAACCTTGTGGCTATGAAGT')),

            Adapter('nCoV-2019_36_LEFT',
                start_sequence=('nCoV-2019_36_LEFT_fw', 'TTAGCTTGGTTGTACGCTGCTG'),
                end_sequence=('nCoV-2019_36_LEFT_revcomp', 'CAGCAGCGTACAACCAAGCTAA')),

            Adapter('nCoV-2019_36_RIGHT',
                start_sequence=('nCoV-2019_36_RIGHT_fw', 'GAACAAAGACCATTGAGTACTCTGGA'),
                end_sequence=('nCoV-2019_36_RIGHT_revcomp', 'TCCAGAGTACTCAATGGTCTTTGTTC')),

            Adapter('nCoV-2019_37_LEFT',
                start_sequence=('nCoV-2019_37_LEFT_fw', 'ACACACCACTGGTTGTTACTCAC'),
                end_sequence=('nCoV-2019_37_LEFT_revcomp', 'GTGAGTAACAACCAGTGGTGTGT')),

            Adapter('nCoV-2019_37_RIGHT',
                start_sequence=('nCoV-2019_37_RIGHT_fw', 'GTCCACACTCTCCTAGCACCAT'),
                end_sequence=('nCoV-2019_37_RIGHT_revcomp', 'ATGGTGCTAGGAGAGTGTGGAC')),

            Adapter('nCoV-2019_38_LEFT',
                start_sequence=('nCoV-2019_38_LEFT_fw', 'ACTGTGTTATGTATGCATCAGCTGT'),
                end_sequence=('nCoV-2019_38_LEFT_revcomp', 'ACAGCTGATGCATACATAACACAGT')),

            Adapter('nCoV-2019_38_RIGHT',
                start_sequence=('nCoV-2019_38_RIGHT_fw', 'CACCAAGAGTCAGTCTAAAGTAGCG'),
                end_sequence=('nCoV-2019_38_RIGHT_revcomp', 'CGCTACTTTAGACTGACTCTTGGTG')),

            Adapter('nCoV-2019_39_LEFT',
                start_sequence=('nCoV-2019_39_LEFT_fw', 'AGTATTGCCCTATTTTCTTCATAACTGGT'),
                end_sequence=('nCoV-2019_39_LEFT_revcomp', 'ACCAGTTATGAAGAAAATAGGGCAATACT')),

            Adapter('nCoV-2019_39_RIGHT',
                start_sequence=('nCoV-2019_39_RIGHT_fw', 'TGTAACTGGACACATTGAGCCC'),
                end_sequence=('nCoV-2019_39_RIGHT_revcomp', 'GGGCTCAATGTGTCCAGTTACA')),

            Adapter('nCoV-2019_40_LEFT',
                start_sequence=('nCoV-2019_40_LEFT_fw', 'TGCACATCAGTAGTCTTACTCTCAGT'),
                end_sequence=('nCoV-2019_40_LEFT_revcomp', 'ACTGAGAGTAAGACTACTGATGTGCA')),

            Adapter('nCoV-2019_40_RIGHT',
                start_sequence=('nCoV-2019_40_RIGHT_fw', 'CATGGCTGCATCACGGTCAAAT'),
                end_sequence=('nCoV-2019_40_RIGHT_revcomp', 'ATTTGACCGTGATGCAGCCATG')),

            Adapter('nCoV-2019_41_LEFT',
                start_sequence=('nCoV-2019_41_LEFT_fw', 'GTTCCCTTCCATCATATGCAGCT'),
                end_sequence=('nCoV-2019_41_LEFT_revcomp', 'AGCTGCATATGATGGAAGGGAAC')),

            Adapter('nCoV-2019_41_RIGHT',
                start_sequence=('nCoV-2019_41_RIGHT_fw', 'TGGTATGACAACCATTAGTTTGGCT'),
                end_sequence=('nCoV-2019_41_RIGHT_revcomp', 'AGCCAAACTAATGGTTGTCATACCA')),

            Adapter('nCoV-2019_42_LEFT',
                start_sequence=('nCoV-2019_42_LEFT_fw', 'TGCAAGAGATGGTTGTGTTCCC'),
                end_sequence=('nCoV-2019_42_LEFT_revcomp', 'GGGAACACAACCATCTCTTGCA')),

            Adapter('nCoV-2019_42_RIGHT',
                start_sequence=('nCoV-2019_42_RIGHT_fw', 'CCTACCTCCCTTTGTTGTGTTGT'),
                end_sequence=('nCoV-2019_42_RIGHT_revcomp', 'ACAACACAACAAAGGGAGGTAGG')),

            Adapter('nCoV-2019_43_LEFT',
                start_sequence=('nCoV-2019_43_LEFT_fw', 'TACGACAGATGTCTTGTGCTGC'),
                end_sequence=('nCoV-2019_43_LEFT_revcomp', 'GCAGCACAAGACATCTGTCGTA')),

            Adapter('nCoV-2019_43_RIGHT',
                start_sequence=('nCoV-2019_43_RIGHT_fw', 'AGCAGCATCTACAGCAAAAGCA'),
                end_sequence=('nCoV-2019_43_RIGHT_revcomp', 'TGCTTTTGCTGTAGATGCTGCT')),

            Adapter('nCoV-2019_44_LEFT',
                start_sequence=('nCoV-2019_44_LEFT_fw', 'TGCCACAGTACGTCTACAAGCT'),
                end_sequence=('nCoV-2019_44_LEFT_revcomp', 'AGCTTGTAGACGTACTGTGGCA')),

            Adapter('nCoV-2019_44_LEFT_alt3',
                start_sequence=('nCoV-2019_44_LEFT_alt3_fw', 'CCACAGTACGTCTACAAGCTGG'),
                end_sequence=('nCoV-2019_44_LEFT_alt3_revcomp', 'CCAGCTTGTAGACGTACTGTGG')),

            Adapter('nCoV-2019_44_RIGHT',
                start_sequence=('nCoV-2019_44_RIGHT_fw', 'AACCTTTCCACATACCGCAGAC'),
                end_sequence=('nCoV-2019_44_RIGHT_revcomp', 'GTCTGCGGTATGTGGAAAGGTT')),

            Adapter('nCoV-2019_44_RIGHT_alt0',
                start_sequence=('nCoV-2019_44_RIGHT_alt0_fw', 'CGCAGACGGTACAGACTGTGTT'),
                end_sequence=('nCoV-2019_44_RIGHT_alt0_revcomp', 'AACACAGTCTGTACCGTCTGCG')),

            Adapter('nCoV-2019_45_LEFT',
                start_sequence=('nCoV-2019_45_LEFT_fw', 'TACCTACAACTTGTGCTAATGACCC'),
                end_sequence=('nCoV-2019_45_LEFT_revcomp', 'GGGTCATTAGCACAAGTTGTAGGTA')),

            Adapter('nCoV-2019_45_LEFT_alt2',
                start_sequence=('nCoV-2019_45_LEFT_alt2_fw', 'AGTATGTACAAATACCTACAACTTGTGCT'),
                end_sequence=('nCoV-2019_45_LEFT_alt2_revcomp', 'AGCACAAGTTGTAGGTATTTGTACATACT')),

            Adapter('nCoV-2019_45_RIGHT',
                start_sequence=('nCoV-2019_45_RIGHT_fw', 'AAATTGTTTCTTCATGTTGGTAGTTAGAGA'),
                end_sequence=('nCoV-2019_45_RIGHT_revcomp', 'TCTCTAACTACCAACATGAAGAAACAATTT')),

            Adapter('nCoV-2019_45_RIGHT_alt7',
                start_sequence=('nCoV-2019_45_RIGHT_alt7_fw', 'TTCATGTTGGTAGTTAGAGAAAGTGTGTC'),
                end_sequence=('nCoV-2019_45_RIGHT_alt7_revcomp', 'GACACACTTTCTCTAACTACCAACATGAA')),

            Adapter('nCoV-2019_46_LEFT',
                start_sequence=('nCoV-2019_46_LEFT_fw', 'TGTCGCTTCCAAGAAAAGGACG'),
                end_sequence=('nCoV-2019_46_LEFT_revcomp', 'CGTCCTTTTCTTGGAAGCGACA')),

            Adapter('nCoV-2019_46_LEFT_alt1',
                start_sequence=('nCoV-2019_46_LEFT_alt1_fw', 'CGCTTCCAAGAAAAGGACGAAGA'),
                end_sequence=('nCoV-2019_46_LEFT_alt1_revcomp', 'TCTTCGTCCTTTTCTTGGAAGCG')),

            Adapter('nCoV-2019_46_RIGHT',
                start_sequence=('nCoV-2019_46_RIGHT_fw', 'CACGTTCACCTAAGTTGGCGTA'),
                end_sequence=('nCoV-2019_46_RIGHT_revcomp', 'TACGCCAACTTAGGTGAACGTG')),

            Adapter('nCoV-2019_46_RIGHT_alt2',
                start_sequence=('nCoV-2019_46_RIGHT_alt2_fw', 'CACGTTCACCTAAGTTGGCGTAT'),
                end_sequence=('nCoV-2019_46_RIGHT_alt2_revcomp', 'ATACGCCAACTTAGGTGAACGTG')),

            Adapter('nCoV-2019_47_LEFT',
                start_sequence=('nCoV-2019_47_LEFT_fw', 'AGGACTGGTATGATTTTGTAGAAAACCC'),
                end_sequence=('nCoV-2019_47_LEFT_revcomp', 'GGGTTTTCTACAAAATCATACCAGTCCT')),

            Adapter('nCoV-2019_47_RIGHT',
                start_sequence=('nCoV-2019_47_RIGHT_fw', 'AATAACGGTCAAAGAGTTTTAACCTCTC'),
                end_sequence=('nCoV-2019_47_RIGHT_revcomp', 'GAGAGGTTAAAACTCTTTGACCGTTATT')),

            Adapter('nCoV-2019_48_LEFT',
                start_sequence=('nCoV-2019_48_LEFT_fw', 'TGTTGACACTGACTTAACAAAGCCT'),
                end_sequence=('nCoV-2019_48_LEFT_revcomp', 'AGGCTTTGTTAAGTCAGTGTCAACA')),

            Adapter('nCoV-2019_48_RIGHT',
                start_sequence=('nCoV-2019_48_RIGHT_fw', 'TAGATTACCAGAAGCAGCGTGC'),
                end_sequence=('nCoV-2019_48_RIGHT_revcomp', 'GCACGCTGCTTCTGGTAATCTA')),

            Adapter('nCoV-2019_49_LEFT',
                start_sequence=('nCoV-2019_49_LEFT_fw', 'AGGAATTACTTGTGTATGCTGCTGA'),
                end_sequence=('nCoV-2019_49_LEFT_revcomp', 'TCAGCAGCATACACAAGTAATTCCT')),

            Adapter('nCoV-2019_49_RIGHT',
                start_sequence=('nCoV-2019_49_RIGHT_fw', 'TGACGATGACTTGGTTAGCATTAATACA'),
                end_sequence=('nCoV-2019_49_RIGHT_revcomp', 'TGTATTAATGCTAACCAAGTCATCGTCA')),

            Adapter('nCoV-2019_50_LEFT',
                start_sequence=('nCoV-2019_50_LEFT_fw', 'GTTGATAAGTACTTTGATTGTTACGATGGT'),
                end_sequence=('nCoV-2019_50_LEFT_revcomp', 'ACCATCGTAACAATCAAAGTACTTATCAAC')),

            Adapter('nCoV-2019_50_RIGHT',
                start_sequence=('nCoV-2019_50_RIGHT_fw', 'TAACATGTTGTGCCAACCACCA'),
                end_sequence=('nCoV-2019_50_RIGHT_revcomp', 'TGGTGGTTGGCACAACATGTTA')),

            Adapter('nCoV-2019_51_LEFT',
                start_sequence=('nCoV-2019_51_LEFT_fw', 'TCAATAGCCGCCACTAGAGGAG'),
                end_sequence=('nCoV-2019_51_LEFT_revcomp', 'CTCCTCTAGTGGCGGCTATTGA')),

            Adapter('nCoV-2019_51_RIGHT',
                start_sequence=('nCoV-2019_51_RIGHT_fw', 'AGTGCATTAACATTGGCCGTGA'),
                end_sequence=('nCoV-2019_51_RIGHT_revcomp', 'TCACGGCCAATGTTAATGCACT')),

            Adapter('nCoV-2019_52_LEFT',
                start_sequence=('nCoV-2019_52_LEFT_fw', 'CATCAGGAGATGCCACAACTGC'),
                end_sequence=('nCoV-2019_52_LEFT_revcomp', 'GCAGTTGTGGCATCTCCTGATG')),

            Adapter('nCoV-2019_52_RIGHT',
                start_sequence=('nCoV-2019_52_RIGHT_fw', 'GTTGAGAGCAAAATTCATGAGGTCC'),
                end_sequence=('nCoV-2019_52_RIGHT_revcomp', 'GGACCTCATGAATTTTGCTCTCAAC')),

            Adapter('nCoV-2019_53_LEFT',
                start_sequence=('nCoV-2019_53_LEFT_fw', 'AGCAAAATGTTGGACTGAGACTGA'),
                end_sequence=('nCoV-2019_53_LEFT_revcomp', 'TCAGTCTCAGTCCAACATTTTGCT')),

            Adapter('nCoV-2019_53_RIGHT',
                start_sequence=('nCoV-2019_53_RIGHT_fw', 'AGCCTCATAAAACTCAGGTTCCC'),
                end_sequence=('nCoV-2019_53_RIGHT_revcomp', 'GGGAACCTGAGTTTTATGAGGCT')),

            Adapter('nCoV-2019_54_LEFT',
                start_sequence=('nCoV-2019_54_LEFT_fw', 'TGAGTTAACAGGACACATGTTAGACA'),
                end_sequence=('nCoV-2019_54_LEFT_revcomp', 'TGTCTAACATGTGTCCTGTTAACTCA')),

            Adapter('nCoV-2019_54_RIGHT',
                start_sequence=('nCoV-2019_54_RIGHT_fw', 'AACCAAAAACTTGTCCATTAGCACA'),
                end_sequence=('nCoV-2019_54_RIGHT_revcomp', 'TGTGCTAATGGACAAGTTTTTGGTT')),

            Adapter('nCoV-2019_55_LEFT',
                start_sequence=('nCoV-2019_55_LEFT_fw', 'ACTCAACTTTACTTAGGAGGTATGAGCT'),
                end_sequence=('nCoV-2019_55_LEFT_revcomp', 'AGCTCATACCTCCTAAGTAAAGTTGAGT')),

            Adapter('nCoV-2019_55_RIGHT',
                start_sequence=('nCoV-2019_55_RIGHT_fw', 'GGTGTACTCTCCTATTTGTACTTTACTGT'),
                end_sequence=('nCoV-2019_55_RIGHT_revcomp', 'ACAGTAAAGTACAAATAGGAGAGTACACC')),

            Adapter('nCoV-2019_56_LEFT',
                start_sequence=('nCoV-2019_56_LEFT_fw', 'ACCTAGACCACCACTTAACCGA'),
                end_sequence=('nCoV-2019_56_LEFT_revcomp', 'TCGGTTAAGTGGTGGTCTAGGT')),

            Adapter('nCoV-2019_56_RIGHT',
                start_sequence=('nCoV-2019_56_RIGHT_fw', 'ACACTATGCGAGCAGAAGGGTA'),
                end_sequence=('nCoV-2019_56_RIGHT_revcomp', 'TACCCTTCTGCTCGCATAGTGT')),

            Adapter('nCoV-2019_57_LEFT',
                start_sequence=('nCoV-2019_57_LEFT_fw', 'ATTCTACACTCCAGGGACCACC'),
                end_sequence=('nCoV-2019_57_LEFT_revcomp', 'GGTGGTCCCTGGAGTGTAGAAT')),

            Adapter('nCoV-2019_57_RIGHT',
                start_sequence=('nCoV-2019_57_RIGHT_fw', 'GTAATTGAGCAGGGTCGCCAAT'),
                end_sequence=('nCoV-2019_57_RIGHT_revcomp', 'ATTGGCGACCCTGCTCAATTAC')),

            Adapter('nCoV-2019_58_LEFT',
                start_sequence=('nCoV-2019_58_LEFT_fw', 'TGATTTGAGTGTTGTCAATGCCAGA'),
                end_sequence=('nCoV-2019_58_LEFT_revcomp', 'TCTGGCATTGACAACACTCAAATCA')),

            Adapter('nCoV-2019_58_RIGHT',
                start_sequence=('nCoV-2019_58_RIGHT_fw', 'CTTTTCTCCAAGCAGGGTTACGT'),
                end_sequence=('nCoV-2019_58_RIGHT_revcomp', 'ACGTAACCCTGCTTGGAGAAAAG')),

            Adapter('nCoV-2019_59_LEFT',
                start_sequence=('nCoV-2019_59_LEFT_fw', 'TCACGCATGATGTTTCATCTGCA'),
                end_sequence=('nCoV-2019_59_LEFT_revcomp', 'TGCAGATGAAACATCATGCGTGA')),

            Adapter('nCoV-2019_59_RIGHT',
                start_sequence=('nCoV-2019_59_RIGHT_fw', 'AAGAGTCCTGTTACATTTTCAGCTTG'),
                end_sequence=('nCoV-2019_59_RIGHT_revcomp', 'CAAGCTGAAAATGTAACAGGACTCTT')),

            Adapter('nCoV-2019_60_LEFT',
                start_sequence=('nCoV-2019_60_LEFT_fw', 'TGATAGAGACCTTTATGACAAGTTGCA'),
                end_sequence=('nCoV-2019_60_LEFT_revcomp', 'TGCAACTTGTCATAAAGGTCTCTATCA')),

            Adapter('nCoV-2019_60_RIGHT',
                start_sequence=('nCoV-2019_60_RIGHT_fw', 'GGTACCAACAGCTTCTCTAGTAGC'),
                end_sequence=('nCoV-2019_60_RIGHT_revcomp', 'GCTACTAGAGAAGCTGTTGGTACC')),

            Adapter('nCoV-2019_61_LEFT',
                start_sequence=('nCoV-2019_61_LEFT_fw', 'TGTTTATCACCCGCGAAGAAGC'),
                end_sequence=('nCoV-2019_61_LEFT_revcomp', 'GCTTCTTCGCGGGTGATAAACA')),

            Adapter('nCoV-2019_61_RIGHT',
                start_sequence=('nCoV-2019_61_RIGHT_fw', 'ATCACATAGACAACAGGTGCGC'),
                end_sequence=('nCoV-2019_61_RIGHT_revcomp', 'GCGCACCTGTTGTCTATGTGAT')),

            Adapter('nCoV-2019_62_LEFT',
                start_sequence=('nCoV-2019_62_LEFT_fw', 'GGCACATGGCTTTGAGTTGACA'),
                end_sequence=('nCoV-2019_62_LEFT_revcomp', 'TGTCAACTCAAAGCCATGTGCC')),

            Adapter('nCoV-2019_62_RIGHT',
                start_sequence=('nCoV-2019_62_RIGHT_fw', 'GTTGAACCTTTCTACAAGCCGC'),
                end_sequence=('nCoV-2019_62_RIGHT_revcomp', 'GCGGCTTGTAGAAAGGTTCAAC')),

            Adapter('nCoV-2019_63_LEFT',
                start_sequence=('nCoV-2019_63_LEFT_fw', 'TGTTAAGCGTGTTGACTGGACT'),
                end_sequence=('nCoV-2019_63_LEFT_revcomp', 'AGTCCAGTCAACACGCTTAACA')),

            Adapter('nCoV-2019_63_RIGHT',
                start_sequence=('nCoV-2019_63_RIGHT_fw', 'ACAAACTGCCACCATCACAACC'),
                end_sequence=('nCoV-2019_63_RIGHT_revcomp', 'GGTTGTGATGGTGGCAGTTTGT')),

            Adapter('nCoV-2019_64_LEFT',
                start_sequence=('nCoV-2019_64_LEFT_fw', 'TCGATAGATATCCTGCTAATTCCATTGT'),
                end_sequence=('nCoV-2019_64_LEFT_revcomp', 'ACAATGGAATTAGCAGGATATCTATCGA')),

            Adapter('nCoV-2019_64_RIGHT',
                start_sequence=('nCoV-2019_64_RIGHT_fw', 'AGTCTTGTAAAAGTGTTCCAGAGGT'),
                end_sequence=('nCoV-2019_64_RIGHT_revcomp', 'ACCTCTGGAACACTTTTACAAGACT')),

            Adapter('nCoV-2019_65_LEFT',
                start_sequence=('nCoV-2019_65_LEFT_fw', 'GCTGGCTTTAGCTTGTGGGTTT'),
                end_sequence=('nCoV-2019_65_LEFT_revcomp', 'AAACCCACAAGCTAAAGCCAGC')),

            Adapter('nCoV-2019_65_RIGHT',
                start_sequence=('nCoV-2019_65_RIGHT_fw', 'TGTCAGTCATAGAACAAACACCAATAGT'),
                end_sequence=('nCoV-2019_65_RIGHT_revcomp', 'ACTATTGGTGTTTGTTCTATGACTGACA')),

            Adapter('nCoV-2019_66_LEFT',
                start_sequence=('nCoV-2019_66_LEFT_fw', 'GGGTGTGGACATTGCTGCTAAT'),
                end_sequence=('nCoV-2019_66_LEFT_revcomp', 'ATTAGCAGCAATGTCCACACCC')),

            Adapter('nCoV-2019_66_RIGHT',
                start_sequence=('nCoV-2019_66_RIGHT_fw', 'TCAATTTCCATTTGACTCCTGGGT'),
                end_sequence=('nCoV-2019_66_RIGHT_revcomp', 'ACCCAGGAGTCAAATGGAAATTGA')),

            Adapter('nCoV-2019_67_LEFT',
                start_sequence=('nCoV-2019_67_LEFT_fw', 'GTTGTCCAACAATTACCTGAAACTTACT'),
                end_sequence=('nCoV-2019_67_LEFT_revcomp', 'AGTAAGTTTCAGGTAATTGTTGGACAAC')),

            Adapter('nCoV-2019_67_RIGHT',
                start_sequence=('nCoV-2019_67_RIGHT_fw', 'CAACCTTAGAAACTACAGATAAATCTTGGG'),
                end_sequence=('nCoV-2019_67_RIGHT_revcomp', 'CCCAAGATTTATCTGTAGTTTCTAAGGTTG')),

            Adapter('nCoV-2019_68_LEFT',
                start_sequence=('nCoV-2019_68_LEFT_fw', 'ACAGGTTCATCTAAGTGTGTGTGT'),
                end_sequence=('nCoV-2019_68_LEFT_revcomp', 'ACACACACACTTAGATGAACCTGT')),

            Adapter('nCoV-2019_68_RIGHT',
                start_sequence=('nCoV-2019_68_RIGHT_fw', 'CTCCTTTATCAGAACCAGCACCA'),
                end_sequence=('nCoV-2019_68_RIGHT_revcomp', 'TGGTGCTGGTTCTGATAAAGGAG')),

            Adapter('nCoV-2019_69_LEFT',
                start_sequence=('nCoV-2019_69_LEFT_fw', 'TGTCGCAAAATATACTCAACTGTGTCA'),
                end_sequence=('nCoV-2019_69_LEFT_revcomp', 'TGACACAGTTGAGTATATTTTGCGACA')),

            Adapter('nCoV-2019_69_RIGHT',
                start_sequence=('nCoV-2019_69_RIGHT_fw', 'TCTTTATAGCCACGGAACCTCCA'),
                end_sequence=('nCoV-2019_69_RIGHT_revcomp', 'TGGAGGTTCCGTGGCTATAAAGA')),

            Adapter('nCoV-2019_70_LEFT',
                start_sequence=('nCoV-2019_70_LEFT_fw', 'ACAAAAGAAAATGACTCTAAAGAGGGTTT'),
                end_sequence=('nCoV-2019_70_LEFT_revcomp', 'AAACCCTCTTTAGAGTCATTTTCTTTTGT')),

            Adapter('nCoV-2019_70_RIGHT',
                start_sequence=('nCoV-2019_70_RIGHT_fw', 'TGACCTTCTTTTAAAGACATAACAGCAG'),
                end_sequence=('nCoV-2019_70_RIGHT_revcomp', 'CTGCTGTTATGTCTTTAAAAGAAGGTCA')),

            Adapter('nCoV-2019_71_LEFT',
                start_sequence=('nCoV-2019_71_LEFT_fw', 'ACAAATCCAATTCAGTTGTCTTCCTATTC'),
                end_sequence=('nCoV-2019_71_LEFT_revcomp', 'GAATAGGAAGACAACTGAATTGGATTTGT')),

            Adapter('nCoV-2019_71_RIGHT',
                start_sequence=('nCoV-2019_71_RIGHT_fw', 'TGGAAAAGAAAGGTAAGAACAAGTCCT'),
                end_sequence=('nCoV-2019_71_RIGHT_revcomp', 'AGGACTTGTTCTTACCTTTCTTTTCCA')),

            Adapter('nCoV-2019_72_LEFT',
                start_sequence=('nCoV-2019_72_LEFT_fw', 'ACACGTGGTGTTTATTACCCTGAC'),
                end_sequence=('nCoV-2019_72_LEFT_revcomp', 'GTCAGGGTAATAAACACCACGTGT')),

            Adapter('nCoV-2019_72_RIGHT',
                start_sequence=('nCoV-2019_72_RIGHT_fw', 'ACTCTGAACTCACTTTCCATCCAAC'),
                end_sequence=('nCoV-2019_72_RIGHT_revcomp', 'GTTGGATGGAAAGTGAGTTCAGAGT')),

            Adapter('nCoV-2019_73_LEFT',
                start_sequence=('nCoV-2019_73_LEFT_fw', 'CAATTTTGTAATGATCCATTTTTGGGTGT'),
                end_sequence=('nCoV-2019_73_LEFT_revcomp', 'ACACCCAAAAATGGATCATTACAAAATTG')),

            Adapter('nCoV-2019_73_RIGHT',
                start_sequence=('nCoV-2019_73_RIGHT_fw', 'CACCAGCTGTCCAACCTGAAGA'),
                end_sequence=('nCoV-2019_73_RIGHT_revcomp', 'TCTTCAGGTTGGACAGCTGGTG')),

            Adapter('nCoV-2019_74_LEFT',
                start_sequence=('nCoV-2019_74_LEFT_fw', 'ACATCACTAGGTTTCAAACTTTACTTGC'),
                end_sequence=('nCoV-2019_74_LEFT_revcomp', 'GCAAGTAAAGTTTGAAACCTAGTGATGT')),

            Adapter('nCoV-2019_74_RIGHT',
                start_sequence=('nCoV-2019_74_RIGHT_fw', 'GCAACACAGTTGCTGATTCTCTTC'),
                end_sequence=('nCoV-2019_74_RIGHT_revcomp', 'GAAGAGAATCAGCAACTGTGTTGC')),

            Adapter('nCoV-2019_75_LEFT',
                start_sequence=('nCoV-2019_75_LEFT_fw', 'AGAGTCCAACCAACAGAATCTATTGT'),
                end_sequence=('nCoV-2019_75_LEFT_revcomp', 'ACAATAGATTCTGTTGGTTGGACTCT')),

            Adapter('nCoV-2019_75_RIGHT',
                start_sequence=('nCoV-2019_75_RIGHT_fw', 'ACCACCAACCTTAGAATCAAGATTGT'),
                end_sequence=('nCoV-2019_75_RIGHT_revcomp', 'ACAATCTTGATTCTAAGGTTGGTGGT')),

            Adapter('nCoV-2019_76_LEFT',
                start_sequence=('nCoV-2019_76_LEFT_fw', 'AGGGCAAACTGGAAAGATTGCT'),
                end_sequence=('nCoV-2019_76_LEFT_revcomp', 'AGCAATCTTTCCAGTTTGCCCT')),

            Adapter('nCoV-2019_76_LEFT_alt3',
                start_sequence=('nCoV-2019_76_LEFT_alt3_fw', 'GGGCAAACTGGAAAGATTGCTGA'),
                end_sequence=('nCoV-2019_76_LEFT_alt3_revcomp', 'TCAGCAATCTTTCCAGTTTGCCC')),

            Adapter('nCoV-2019_76_RIGHT',
                start_sequence=('nCoV-2019_76_RIGHT_fw', 'ACACCTGTGCCTGTTAAACCAT'),
                end_sequence=('nCoV-2019_76_RIGHT_revcomp', 'ATGGTTTAACAGGCACAGGTGT')),

            Adapter('nCoV-2019_76_RIGHT_alt0',
                start_sequence=('nCoV-2019_76_RIGHT_alt0_fw', 'ACCTGTGCCTGTTAAACCATTGA'),
                end_sequence=('nCoV-2019_76_RIGHT_alt0_revcomp', 'TCAATGGTTTAACAGGCACAGGT')),

            Adapter('nCoV-2019_77_LEFT',
                start_sequence=('nCoV-2019_77_LEFT_fw', 'CCAGCAACTGTTTGTGGACCTA'),
                end_sequence=('nCoV-2019_77_LEFT_revcomp', 'TAGGTCCACAAACAGTTGCTGG')),

            Adapter('nCoV-2019_77_RIGHT',
                start_sequence=('nCoV-2019_77_RIGHT_fw', 'CAGCCCCTATTAAACAGCCTGC'),
                end_sequence=('nCoV-2019_77_RIGHT_revcomp', 'GCAGGCTGTTTAATAGGGGCTG')),

            Adapter('nCoV-2019_78_LEFT',
                start_sequence=('nCoV-2019_78_LEFT_fw', 'CAACTTACTCCTACTTGGCGTGT'),
                end_sequence=('nCoV-2019_78_LEFT_revcomp', 'ACACGCCAAGTAGGAGTAAGTTG')),

            Adapter('nCoV-2019_78_RIGHT',
                start_sequence=('nCoV-2019_78_RIGHT_fw', 'TGTGTACAAAAACTGCCATATTGCA'),
                end_sequence=('nCoV-2019_78_RIGHT_revcomp', 'TGCAATATGGCAGTTTTTGTACACA')),

            Adapter('nCoV-2019_79_LEFT',
                start_sequence=('nCoV-2019_79_LEFT_fw', 'GTGGTGATTCAACTGAATGCAGC'),
                end_sequence=('nCoV-2019_79_LEFT_revcomp', 'GCTGCATTCAGTTGAATCACCAC')),

            Adapter('nCoV-2019_79_RIGHT',
                start_sequence=('nCoV-2019_79_RIGHT_fw', 'CATTTCATCTGTGAGCAAAGGTGG'),
                end_sequence=('nCoV-2019_79_RIGHT_revcomp', 'CCACCTTTGCTCACAGATGAAATG')),

            Adapter('nCoV-2019_80_LEFT',
                start_sequence=('nCoV-2019_80_LEFT_fw', 'TTGCCTTGGTGATATTGCTGCT'),
                end_sequence=('nCoV-2019_80_LEFT_revcomp', 'AGCAGCAATATCACCAAGGCAA')),

            Adapter('nCoV-2019_80_RIGHT',
                start_sequence=('nCoV-2019_80_RIGHT_fw', 'TGGAGCTAAGTTGTTTAACAAGCG'),
                end_sequence=('nCoV-2019_80_RIGHT_revcomp', 'CGCTTGTTAAACAACTTAGCTCCA')),

            Adapter('nCoV-2019_81_LEFT',
                start_sequence=('nCoV-2019_81_LEFT_fw', 'GCACTTGGAAAACTTCAAGATGTGG'),
                end_sequence=('nCoV-2019_81_LEFT_revcomp', 'CCACATCTTGAAGTTTTCCAAGTGC')),

            Adapter('nCoV-2019_81_RIGHT',
                start_sequence=('nCoV-2019_81_RIGHT_fw', 'GTGAAGTTCTTTTCTTGTGCAGGG'),
                end_sequence=('nCoV-2019_81_RIGHT_revcomp', 'CCCTGCACAAGAAAAGAACTTCAC')),

            Adapter('nCoV-2019_82_LEFT',
                start_sequence=('nCoV-2019_82_LEFT_fw', 'GGGCTATCATCTTATGTCCTTCCCT'),
                end_sequence=('nCoV-2019_82_LEFT_revcomp', 'AGGGAAGGACATAAGATGATAGCCC')),

            Adapter('nCoV-2019_82_RIGHT',
                start_sequence=('nCoV-2019_82_RIGHT_fw', 'TGCCAGAGATGTCACCTAAATCAA'),
                end_sequence=('nCoV-2019_82_RIGHT_revcomp', 'TTGATTTAGGTGACATCTCTGGCA')),

            Adapter('nCoV-2019_83_LEFT',
                start_sequence=('nCoV-2019_83_LEFT_fw', 'TCCTTTGCAACCTGAATTAGACTCA'),
                end_sequence=('nCoV-2019_83_LEFT_revcomp', 'TGAGTCTAATTCAGGTTGCAAAGGA')),

            Adapter('nCoV-2019_83_RIGHT',
                start_sequence=('nCoV-2019_83_RIGHT_fw', 'TTTGACTCCTTTGAGCACTGGC'),
                end_sequence=('nCoV-2019_83_RIGHT_revcomp', 'GCCAGTGCTCAAAGGAGTCAAA')),

            Adapter('nCoV-2019_84_LEFT',
                start_sequence=('nCoV-2019_84_LEFT_fw', 'TGCTGTAGTTGTCTCAAGGGCT'),
                end_sequence=('nCoV-2019_84_LEFT_revcomp', 'AGCCCTTGAGACAACTACAGCA')),

            Adapter('nCoV-2019_84_RIGHT',
                start_sequence=('nCoV-2019_84_RIGHT_fw', 'AGGTGTGAGTAAACTGTTACAAACAAC'),
                end_sequence=('nCoV-2019_84_RIGHT_revcomp', 'GTTGTTTGTAACAGTTTACTCACACCT')),

            Adapter('nCoV-2019_85_LEFT',
                start_sequence=('nCoV-2019_85_LEFT_fw', 'ACTAGCACTCTCCAAGGGTGTT'),
                end_sequence=('nCoV-2019_85_LEFT_revcomp', 'AACACCCTTGGAGAGTGCTAGT')),

            Adapter('nCoV-2019_85_RIGHT',
                start_sequence=('nCoV-2019_85_RIGHT_fw', 'ACACAGTCTTTTACTCCAGATTCCC'),
                end_sequence=('nCoV-2019_85_RIGHT_revcomp', 'GGGAATCTGGAGTAAAAGACTGTGT')),

            Adapter('nCoV-2019_86_LEFT',
                start_sequence=('nCoV-2019_86_LEFT_fw', 'TCAGGTGATGGCACAACAAGTC'),
                end_sequence=('nCoV-2019_86_LEFT_revcomp', 'GACTTGTTGTGCCATCACCTGA')),

            Adapter('nCoV-2019_86_RIGHT',
                start_sequence=('nCoV-2019_86_RIGHT_fw', 'ACGAAAGCAAGAAAAAGAAGTACGC'),
                end_sequence=('nCoV-2019_86_RIGHT_revcomp', 'GCGTACTTCTTTTTCTTGCTTTCGT')),

            Adapter('nCoV-2019_87_LEFT',
                start_sequence=('nCoV-2019_87_LEFT_fw', 'CGACTACTAGCGTGCCTTTGTA'),
                end_sequence=('nCoV-2019_87_LEFT_revcomp', 'TACAAAGGCACGCTAGTAGTCG')),

            Adapter('nCoV-2019_87_RIGHT',
                start_sequence=('nCoV-2019_87_RIGHT_fw', 'ACTAGGTTCCATTGTTCAAGGAGC'),
                end_sequence=('nCoV-2019_87_RIGHT_revcomp', 'GCTCCTTGAACAATGGAACCTAGT')),

            Adapter('nCoV-2019_88_LEFT',
                start_sequence=('nCoV-2019_88_LEFT_fw', 'CCATGGCAGATTCCAACGGTAC'),
                end_sequence=('nCoV-2019_88_LEFT_revcomp', 'GTACCGTTGGAATCTGCCATGG')),

            Adapter('nCoV-2019_88_RIGHT',
                start_sequence=('nCoV-2019_88_RIGHT_fw', 'TGGTCAGAATAGTGCCATGGAGT'),
                end_sequence=('nCoV-2019_88_RIGHT_revcomp', 'ACTCCATGGCACTATTCTGACCA')),

            Adapter('nCoV-2019_89_LEFT',
                start_sequence=('nCoV-2019_89_LEFT_fw', 'GTACGCGTTCCATGTGGTCATT'),
                end_sequence=('nCoV-2019_89_LEFT_revcomp', 'AATGACCACATGGAACGCGTAC')),

            Adapter('nCoV-2019_89_LEFT_alt2',
                start_sequence=('nCoV-2019_89_LEFT_alt2_fw', 'CGCGTTCCATGTGGTCATTCAA'),
                end_sequence=('nCoV-2019_89_LEFT_alt2_revcomp', 'TTGAATGACCACATGGAACGCG')),

            Adapter('nCoV-2019_89_RIGHT',
                start_sequence=('nCoV-2019_89_RIGHT_fw', 'ACCTGAAAGTCAACGAGATGAAACA'),
                end_sequence=('nCoV-2019_89_RIGHT_revcomp', 'TGTTTCATCTCGTTGACTTTCAGGT')),

            Adapter('nCoV-2019_89_RIGHT_alt4',
                start_sequence=('nCoV-2019_89_RIGHT_alt4_fw', 'ACGAGATGAAACATCTGTTGTCACT'),
                end_sequence=('nCoV-2019_89_RIGHT_alt4_revcomp', 'AGTGACAACAGATGTTTCATCTCGT')),

            Adapter('nCoV-2019_90_LEFT',
                start_sequence=('nCoV-2019_90_LEFT_fw', 'ACACAGACCATTCCAGTAGCAGT'),
                end_sequence=('nCoV-2019_90_LEFT_revcomp', 'ACTGCTACTGGAATGGTCTGTGT')),

            Adapter('nCoV-2019_90_RIGHT',
                start_sequence=('nCoV-2019_90_RIGHT_fw', 'TGAAATGGTGAATTGCCCTCGT'),
                end_sequence=('nCoV-2019_90_RIGHT_revcomp', 'ACGAGGGCAATTCACCATTTCA')),

            Adapter('nCoV-2019_91_LEFT',
                start_sequence=('nCoV-2019_91_LEFT_fw', 'TCACTACCAAGAGTGTGTTAGAGGT'),
                end_sequence=('nCoV-2019_91_LEFT_revcomp', 'ACCTCTAACACACTCTTGGTAGTGA')),

            Adapter('nCoV-2019_91_RIGHT',
                start_sequence=('nCoV-2019_91_RIGHT_fw', 'TTCAAGTGAGAACCAAAAGATAATAAGCA'),
                end_sequence=('nCoV-2019_91_RIGHT_revcomp', 'TGCTTATTATCTTTTGGTTCTCACTTGAA')),

            Adapter('nCoV-2019_92_LEFT',
                start_sequence=('nCoV-2019_92_LEFT_fw', 'TTTGTGCTTTTTAGCCTTTCTGCT'),
                end_sequence=('nCoV-2019_92_LEFT_revcomp', 'AGCAGAAAGGCTAAAAAGCACAAA')),

            Adapter('nCoV-2019_92_RIGHT',
                start_sequence=('nCoV-2019_92_RIGHT_fw', 'AGGTTCCTGGCAATTAATTGTAAAAGG'),
                end_sequence=('nCoV-2019_92_RIGHT_revcomp', 'CCTTTTACAATTAATTGCCAGGAACCT')),

            Adapter('nCoV-2019_93_LEFT',
                start_sequence=('nCoV-2019_93_LEFT_fw', 'TGAGGCTGGTTCTAAATCACCCA'),
                end_sequence=('nCoV-2019_93_LEFT_revcomp', 'TGGGTGATTTAGAACCAGCCTCA')),

            Adapter('nCoV-2019_93_RIGHT',
                start_sequence=('nCoV-2019_93_RIGHT_fw', 'AGGTCTTCCTTGCCATGTTGAG'),
                end_sequence=('nCoV-2019_93_RIGHT_revcomp', 'CTCAACATGGCAAGGAAGACCT')),

            Adapter('nCoV-2019_94_LEFT',
                start_sequence=('nCoV-2019_94_LEFT_fw', 'GGCCCCAAGGTTTACCCAATAA'),
                end_sequence=('nCoV-2019_94_LEFT_revcomp', 'TTATTGGGTAAACCTTGGGGCC')),

            Adapter('nCoV-2019_94_RIGHT',
                start_sequence=('nCoV-2019_94_RIGHT_fw', 'TTTGGCAATGTTGTTCCTTGAGG'),
                end_sequence=('nCoV-2019_94_RIGHT_revcomp', 'CCTCAAGGAACAACATTGCCAAA')),

            Adapter('nCoV-2019_95_LEFT',
                start_sequence=('nCoV-2019_95_LEFT_fw', 'TGAGGGAGCCTTGAATACACCA'),
                end_sequence=('nCoV-2019_95_LEFT_revcomp', 'TGGTGTATTCAAGGCTCCCTCA')),

            Adapter('nCoV-2019_95_RIGHT',
                start_sequence=('nCoV-2019_95_RIGHT_fw', 'CAGTACGTTTTTGCCGAGGCTT'),
                end_sequence=('nCoV-2019_95_RIGHT_revcomp', 'AAGCCTCGGCAAAAACGTACTG')),

            Adapter('nCoV-2019_96_LEFT',
                start_sequence=('nCoV-2019_96_LEFT_fw', 'GCCAACAACAACAAGGCCAAAC'),
                end_sequence=('nCoV-2019_96_LEFT_revcomp', 'GTTTGGCCTTGTTGTTGTTGGC')),

            Adapter('nCoV-2019_96_RIGHT',
                start_sequence=('nCoV-2019_96_RIGHT_fw', 'TAGGCTCTGTTGGTGGGAATGT'),
                end_sequence=('nCoV-2019_96_RIGHT_revcomp', 'ACATTCCCACCAACAGAGCCTA')),

            Adapter('nCoV-2019_97_LEFT',
                start_sequence=('nCoV-2019_97_LEFT_fw', 'TGGATGACAAAGATCCAAATTTCAAAGA'),
                end_sequence=('nCoV-2019_97_LEFT_revcomp', 'TCTTTGAAATTTGGATCTTTGTCATCCA')),

            Adapter('nCoV-2019_97_RIGHT',
                start_sequence=('nCoV-2019_97_RIGHT_fw', 'ACACACTGATTAAAGATTGCTATGTGAG'),
                end_sequence=('nCoV-2019_97_RIGHT_revcomp', 'CTCACATAGCAATCTTTAATCAGTGTGT')),

            Adapter('nCoV-2019_98_LEFT',
                start_sequence=('nCoV-2019_98_LEFT_fw', 'AACAATTGCAACAATCCATGAGCA'),
                end_sequence=('nCoV-2019_98_LEFT_revcomp', 'TGCTCATGGATTGTTGCAATTGTT')),

            Adapter('nCoV-2019_98_RIGHT',
                start_sequence=('nCoV-2019_98_RIGHT_fw', 'TTCTCCTAAGAAGCTATTAAAATCACATGG'),
                end_sequence=('nCoV-2019_98_RIGHT_revcomp', 'CCATGTGATTTTAATAGCTTCTTAGGAGAA')),
            ]


def make_full_native_barcode_adapter(barcode_num):
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (reverse)'][0]
    start_barcode_seq = barcode.start_sequence[1]
    end_barcode_seq = barcode.end_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA' + start_barcode_seq + 'CAGCACCT'
    end_full_seq = 'AGGTGCTG' + end_barcode_seq + 'TTAACCTTAGCAATACGTAACTGAACGAAGT'

    return Adapter('Native barcoding ' + str(barcode_num) + ' (full sequence)',
                   start_sequence=('NB' + '%02d' % barcode_num + '_start', start_full_seq),
                   end_sequence=('NB' + '%02d' % barcode_num + '_end', end_full_seq))


def make_old_full_rapid_barcode_adapter(barcode_num):  # applies to SQK-RBK001
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACG' + 'TATTGCT' + start_barcode_seq + \
                     'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence, old)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))


def make_new_full_rapid_barcode_adapter(barcode_num):  # applies to SQK-RBK004
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACG' + 'GCTTGGGTGTTTAACC' + start_barcode_seq + \
                     'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence, new)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))
