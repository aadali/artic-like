from sys import argv

from functions import call_snp

input_bam, genome, output = argv[1:]

if __name__ == '__main__':
    input_bam, genome, output = argv[1:]
    call_snp(input_bam, genome, output)
