from sys import argv
from os import path
import argparse
import vcf
import re

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


pass_vcf = open(passvcf, "w")
fail_vcf = open(failvcf, "w")
with open(inputvcf, "r") as invcf:
    for line in invcf:
        if line.startswith("#"):
            pass_vcf.write(line)
            fail_vcf.write(line)
            continue
        CHRO, POS, _, REF, ALT, QUAL, _, INFO, FORMAT, SAMPLE = line.strip().split("\t")
        if len(ALT.split(",")) > 1:
            raise Exception("This code doesn't support multiple genotypes")
        bases = len(REF) - len(ALT)
        if bases % 3 == 0:
            inframe = True
        else:
            inframe = False
        GT, GQ = SAMPLE.split(":")
        GT = int(GT)
        GQ = int(GQ)
        dp = int(re.search("DP=(\d+)", INFO).group(1))
        if GT != 1:
            raise Exception("Medaka variants is not alt homo")
        if GQ < args.min_qual or float(QUAL) < args.min_qual or int(dp) < args.min_depth:
            fail_vcf.write(line)
        else:
            pass_vcf.write(line)

        continue # don't run the below code

        gt_idx = FORMAT.split(":").index("GT")
        sample = SAMPLE.split(":")
        gt = sample[gt_idx]
        dp = int(re.search("DP=(\d+)", INFO).group(1))

        # the low dp or low QUAL record will be treated as bad variant
        is_failed = dp < args.min_depth or float(QUAL) < args.min_qual or (not args.frameshifts and not inframe)

        if not (dp < args.min_depth or float(QUAL) < args.min_qual):
            # reportvcf contains the good record ignore its heterozygosity, will be displayed in the pdf report
            report_vcf.write(line)
        if medaka_vcf:
            is_het = False if gt == "1" else True
            # the variants in medaka vcf must be ref_hom or alt_hom. het genotype is impossible
            if is_failed:
                fail_vcf.write(line)
            else:
                pass_vcf.write(line)

        else:
            # for longshot vcf
            is_het = True if gt == "0/1" else False
            if is_het:
                if args.het_site == "more":
                    ac = re.search("AC=(\d+),(\d+)", INFO)
                    ref_ac = int(ac.group(1))
                    alt_ac = int(ac.group(2))
                    if is_failed:
                        fail_vcf.write(line)
                    else:
                        # if this variant is a good het site and the alt dp is bigger than ref dp, then this record will
                        # be write into the pass vcf, and the pass vcf will be called by bcftools consensus. And bcftools
                        # consensus will apply the alt allele ignore its genotype. this is what "het_site: more" means
                        # But if alt dp is less than ref_ac, then ignore this site, and the position of final consensus will
                        # be ref allele. Although this position is ignored, it can be found in the report vcf
                        if alt_ac > ref_ac:
                            pass_vcf.write(line)

                # if specify the het_site in config.yaml == "N", all the het sites will be N in the final consensus
                elif args.het_site == "N":
                    fail_vcf.write(line)
            else:
                if is_failed:
                    fail_vcf.write(line)
                else:
                    pass_vcf.write(line)
