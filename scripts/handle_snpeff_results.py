from sys import argv
from os import path
from functions import write_annotate_result

usage = f"usage: {path.basename(__file__)} <invcf> <outfile>"

if __name__ == '__main__':
    if len(argv) != 3:
        raise Exception(usage)
    else:
        write_annotate_result(argv[1], argv[2])
