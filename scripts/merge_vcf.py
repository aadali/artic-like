from sys import argv
from os import path

from functions import merge

if __name__ == '__main__':

    usage = f"python3 {path.basename(__file__)} <counts> <genome.fa> " \
            f"<input_longhshot_ann_vcf> <output_longshot_vcf> " \
            f"<input_medaka_pass_vcf> <output_medaka_pass_vcf>"

    if len(argv) != 5:
        raise Exception(usage)
    inputvcf, counts, genome, outvcf = argv[1:]
    merge(inputvcf, counts, genome, outvcf)
