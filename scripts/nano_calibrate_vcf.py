import pandas as pd
from sys import argv
from os import path
import gzip

"""
select the snv whose freq is great than min_freq as variant from the results of igvtools counts by window=1.
if the ref_ac > alt_ac, the position in consensus will be ref.
elif the ref_ac < alt_ac, the position in consensus will be alt.
for deletion range, the script will ignore it.
"""

min_freq = 0.1

usage = f"python {path.basename(__file__)} <igv_counts> <genome.fa> " \
        f"<longhshot_ann_vcf> <modify_longshot_out_vcf> " \
        f"<medaka_pass_vcf> <medaka_pass_out_vcf>"

if len(argv) != 7:
    raise Exception(usage)
# igv_counts = "/home/a/big/ycq/ncov_results/anyang_bar03/anyang_bar03/aligns/out.txt"
# longshot_ann = "/home/a/big/ycq/ncov_results/anyang_bar03/anyang_bar03/variants/anyang_bar03.longshot.ann.vcf"
# medaka_pass = "/home/a/big/ycq/ncov_results/anyang_bar03/anyang_bar03/variants/anyang_bar03.pass.vcf"
# ref_genome = "/home/a/big/ycq/projects/artic-like/genome/sars-cov-2/sequences.fa"
igv_counts, ref_genome, longshot_ann, longshot_ann_out, medaka_pass, medaka_pass_out = argv[1:]


def func(vcffile):
    """
    get the vcf's header, dels and the del range
    :param vcffile:
    :return:
    """
    headers = []
    dels = []
    del_range = []
    with open(vcffile, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                headers.append(line.strip())
            else:
                fileds = line.strip().split("\t")
                ref, alt = fileds[3:5]
                pos = int(fileds[1])
                if len(ref) <= len(alt):
                    continue
                else:
                    del_range.extend(list(range(pos, pos + len(ref))))
                    dels.append(fileds)
    return headers, dels, del_range


# longshot vcf is used for reprot dislpaying
longshot_header, longshot_dels, longshot_del_range = func(longshot_ann)

# medaka pass vcf is used for bcftools to generate consensus fasta
medaka_pass_header, medaka_pass_dels, medaka_pass_del_range = func(medaka_pass)
with open(ref_genome, 'r') as genome:
    contig_name = genome.readline().strip()[1:]
    sequences = genome.read().replace("\n", "")

ref_d = dict(zip(range(1, len(sequences) + 1), list(sequences)))

counts_df = pd.read_csv(igv_counts, sep="\t", header=None, usecols=range(5), skiprows=5,
                        names=['pos', 'A', 'C', 'G', 'T'],
                        dtype={'pos': int, 'A': int, 'C': int, 'T': int})

candidate_variants = []
for row in counts_df.itertuples(index=False):
    pos = row.pos
    ref = ref_d[pos]
    base_counts = dict(zip(list("ACGT"), row[1:]))
    sorted_base_counts = list(sorted(base_counts.items(), key=lambda x: x[1], reverse=True))
    total = sum(row[1:]) if sum(row[1:]) else 1
    if ref == sorted_base_counts[0][0]:
        ref, ref_ac = sorted_base_counts[0]
        alt, alt_ac = sorted_base_counts[1]
    elif ref == sorted_base_counts[1][0]:
        ref, ref_ac = sorted_base_counts[1]
        alt, alt_ac = sorted_base_counts[0]
    elif ref == sorted_base_counts[2][0]:
        ref, ref_ac = sorted_base_counts[2]
        alt, alt_ac = sorted_base_counts[0]
    else:
        ref, ref_ac = sorted_base_counts[3]
        alt, alt_ac = sorted_base_counts[0]

    if alt_ac / total > min_freq:
        gt = "1/1" if alt_ac > ref_ac else "0/1"
        # only output the positon which alt_ac >= ref_ac, bcftools will apply this positions to get consensus
        if alt_ac / (alt_ac + ref_ac) > 0.55 and pos not in medaka_pass_del_range:
            var_record = [contig_name, str(pos), ".", ref, alt, "500", 'PASS',
                          f"DP={total};AC={int(ref_ac)},{int(alt_ac)};",
                          "GT:GQ",
                          f"{gt}:500"]
            candidate_variants.append(var_record)

medaka_pass_all = candidate_variants + medaka_pass_dels
longshot_all = candidate_variants + longshot_dels
medaka_pass_all = list(sorted(medaka_pass_all, key=lambda x: int(x[1])))
lonshot_all = list(sorted(longshot_all, key=lambda x: int(x[1])))
with open(longshot_ann_out, 'w') as longshot_out, open(medaka_pass_out, 'w') as medaka_pass_out:
    longshot_out.write("\n".join(longshot_header) + "\n")
    medaka_pass_out.write("\n".join(medaka_pass_header) + "\n")
    longshot_all = list(map(lambda x: "\t".join(x), lonshot_all))
    medaka_pass_all = list(map(lambda x: "\t".join(x), medaka_pass_all))
    longshot_out.write("\n".join(longshot_all))
    medaka_pass_out.write("\n".join(medaka_pass_all))
