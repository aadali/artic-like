from sys import argv
from os import path


def parse_reference(ref_fh):
    contig_seq = dict()
    with open(ref_fh, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                contig = line.strip().split(" ")[0][1:]
                contig_seq[contig] = []
            else:
                contig_seq[contig].append(line)

    for key, value in contig_seq.items():
        contig_seq[key] = list("".join(value).replace("\n", ""))

    return contig_seq


def pre_mask_reference(ref_fh, low_coverage_bed, vcf, out_consensus):
    import gzip
    contig_seq = parse_reference(ref_fh)
    with open(low_coverage_bed, 'r') as low_cov:
        for line in low_cov:
            contig, start, end = line.strip().split("\t")
            for i in range(int(start), int(end)):
                contig_seq[contig][i] = "N"

    for line in gzip.open(vcf, 'rt'):
        if line.startswith("#"):
            continue

        chro, pos, _, ref, alt = line.split("\t")[:5]

        for i in range(0, len(ref)):
            contig_seq[chro][int(pos) - 1 + i] = "N"

    with open(out_consensus, "w") as outfile:
        for contig, seq in contig_seq.items():
            outfile.write(f">{contig}\n")
            outfile.write("".join(seq) + "\n")


if __name__ == '__main__':
    usage = f"usage: {path.basename(__file__)} <reference.fa> <low_coverage.bed> <fail.vcf> <output.consensus>"
    if len(argv) != 5:
        raise Exception(usage + "\n")
    pre_mask_reference(*(argv[1:]))

