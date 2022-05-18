import argparse
from functions import filt_vcf

parser = argparse.ArgumentParser()
parser.add_argument("--min-depth", required=False, type=int, default=20)
parser.add_argument("--min-qual", required=False, type=int, default=20)
parser.add_argument("inputvcf")
parser.add_argument("passvcf")
parser.add_argument("failvcf")

args = parser.parse_args()

inputvcf = args.inputvcf
passvcf = args.passvcf
failvcf = args.failvcf
min_dp = args.min_depth
min_qual = args.min_qual
filt_vcf(inputvcf, passvcf, failvcf, min_dp, min_qual)
