from sys import argv
from os import path

from functions import trim_bam

usage = f"usage: {path.basename(__file__)} <input_bam> <output_bam> <primer_bed> <trim_log> <min_overlap> <trim_bases>"

if __name__ == '__main__':
    paras_num = len(argv)
    if paras_num != 7:
        raise Exception(usage)
    else:
        trim_bam(*argv[1:])
