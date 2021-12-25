from sys import argv
from os import path

import vcf

usage = f"usage: {path.basename(__file__)} <in.vcf> <min_dp> <pass.vcf> <fail.vcf>"

if len(argv) != 5:
    raise Exception(usage + "\n")

vcf_fh = argv[1]
min_dp = argv[2]
pass_vcf_fh = argv[3]
fail_vcf_fh = argv[4]
# vcf_fh = "/home/a/big/ycq/projects/artic-like/test001/variants/test001.longshot.vcf"
invcf = vcf.Reader(filename=vcf_fh)
invcf.filters['low_qual'] = vcf.parser._Filter(id="low_qual", desc="Qual less than 20")
invcf.filters['low_dp'] = vcf.parser._Filter(id="low_dp", desc=f"DP less than {min_dp}")
invcf.filters['het'] = vcf.parser._Filter(id="het", desc="this position is Het")
invcf.filters['fs'] = vcf.parser._Filter(id="fs", desc="frameshift")
pass_vcf = vcf.Writer(open(pass_vcf_fh, "w"), template=invcf)
fail_vcf = vcf.Writer(open(fail_vcf_fh, "w"), template=invcf)


def in_frame(v):
    if len(v.ALT) > 1:
        raise Exception("This code doesn't support multiple genotypes")
    ref, alt = v.REF, v.ALT[0]
    bases = len(ref) - len(alt)
    if bases % 3 == 0:
        return True
    return False


for v in invcf:
    dp = v.INFO['DP']
    qual = v.QUAL

    if not in_frame(v):
        v.FILTER = "fs"
        fail_vcf.write_record(v)
    elif dp < int(min_dp):
        v.FILTER = "low_dp"
        fail_vcf.write_record(v)
    elif qual < 20:
        v.FILTER = "low_qual"
        fail_vcf.write_record(v)
    elif v.num_het:
        v.FILTER = "het"
        fail_vcf.write_record(v)
    else:
        pass_vcf.write_record(v)
