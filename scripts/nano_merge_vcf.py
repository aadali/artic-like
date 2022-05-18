from sys import argv
from os import path

from functions import merge

if __name__ == '__main__':

    usage = f"python3 {path.basename(__file__)} <counts> <genome.fa> " \
            f"<input_longhshot_ann_vcf> <output_longshot_vcf> " \
            f"<input_medaka_pass_vcf> <output_medaka_pass_vcf>"

    if len(argv) != 7:
        raise Exception(usage)
    counts, ref_genome, longshot_ann, longshot_ann_out, medaka_pass, medaka_pass_out = argv[1:]
    merge(longshot_ann, counts, ref_genome, longshot_ann_out)
    merge(medaka_pass, counts, ref_genome, medaka_pass_out)
