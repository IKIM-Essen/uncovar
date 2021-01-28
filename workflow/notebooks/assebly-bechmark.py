import pandas as pd
import glob
import pysam

snakemake_input = glob.glob("../../results/benchmarking/assembly/*")
# snakemake_output = '../../conda'
snakemake_output = "/Users/alex/COVID_Projects/snakemake-workflow-sars-cov2/results/results.csv"

with open(snakemake_output, "w") as out:
# with open(snakemake.output[0], "w") as out:
   print("Accession", "Contigs", 'Total contigs', 'Contig length', 'Reference length', 'Contig frac', 'Edit dist.', 'Edit frac', sep="\t")
   print("Accession", "Contigs", 'Total contigs', 'Contig length', 'Reference length', 'Contig frac', 'Edit dist.', 'Edit frac', sep="\t", file=out)
   sum_of_edit_dist = 0

   for bam_file in snakemake_input:
   # for bam_file in snakemake.input:
      current_contig = 1

      with pysam.AlignmentFile(bam_file, "rb") as samfile:
         total_contigs = samfile.count()
         accession = samfile.get_reference_name(0)

      with pysam.AlignmentFile(bam_file, "rb") as samfile:
         ref_lengths = samfile.lengths[0]
         for read in samfile.fetch():
            query_alignment_length = read.query_alignment_length
            frac = round(query_alignment_length / ref_lengths, 2)
            edit = read.get_tag("NM")
            sum_of_edit_dist = sum_of_edit_dist + int(edit)
            edit_frac = round(read.get_tag("NM")/query_alignment_length, 5)

            print(accession, current_contig, total_contigs, query_alignment_length, ref_lengths, frac, edit, edit_frac, sep="\t")
            print(accession, current_contig, total_contigs, query_alignment_length, ref_lengths, frac, edit, edit_frac, sep="\t", file=out)

            current_contig += 1
   
   print(sum_of_edit_dist)
   print(sum_of_edit_dist, file=out)    