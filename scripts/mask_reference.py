import sys
from sys import argv
from os import path
from Bio import SeqIO
import pandas as pd


def parse_reference(ref_fh):
    with open(ref_fh, 'r') as reference:
        contigs = SeqIO.parse(reference, "fasta")
        return {contig.name: list(contig.seq) for contig in contigs}


def pre_mask_reference(ref_fh, low_coverage_bed, per_pos_each_base_dp, out_consensus):
    contig_seq = parse_reference(ref_fh)
    counts_df = pd.read_csv(per_pos_each_base_dp,
                            sep="\t", usecols=range(6), comment="#", header=None,
                            names=['contig', 'pos', 'A', 'C', 'G', 'T'],
                            dtype={'contig': str, 'pos': int, 'A': int, 'C': int, 'G': int, 'T': int})

    # find the N in reference
    n_idx_in_ref = {}
    for contig, seq in contig_seq.items():
        s = pd.Series(seq)
        ns = s[(s == 'n') | (s == "N")]
        if not ns.empty:
            n_idx_in_ref[contig] = list(ns.index)

    # change the Ns in reference to real bases depend igvtools
    for contig, poses in n_idx_in_ref.items():
        for idx in poses:
            sub_df = counts_df.query("contig == @contig and pos == @idx+1").reset_index(drop=True)
            if sub_df.empty:
                print(f"{contig} at postion {idx} can't find", file=sys.stderr)
                exit(2)
            row = next(sub_df.itertuples(index=False))
            bases_counts = dict(zip(list("ACGT"), row[2:]))
            sorted_base_counts = list(sorted(bases_counts.items(), key=lambda x:x[1], reverse=True))
            contig_seq[contig][idx] = sorted_base_counts[0][0]

    # mask low coverage region
    with open(low_coverage_bed, 'r') as low_cov:
        for line in low_cov:
            contig, start, end = line.strip().split("\t")
            for i in range(int(start), int(end)):
                contig_seq[contig][i] = "N"

    # write out masked reference
    with open(out_consensus, "w") as outfile:
        for contig, seq in contig_seq.items():
            outfile.write(f">{contig}\n")
            outfile.write("".join(seq) + "\n")


if __name__ == '__main__':
    usage = f"usage: {path.basename(__file__)} <reference.fa> <low_coverage.bed> <per_pos_each_base_dp> <output.consensus>"
    if len(argv) != 5:
        raise Exception(usage + "\n")
    pre_mask_reference(*(argv[1:]))
