import vcf
import re
from sys import argv
from os import path

usage = f"usage: {path.basename(__file__)} <invcf> <outfile>"


def write_annotate_result(invcf, outfile):
    with open(outfile, "w") as outf:
        variants = vcf.Reader(open(invcf, "r"))
        mat = re.search("'(.*)'", variants.infos['ANN'].desc)
        ann_keys = mat.group(1).split("|")
        ann_keys = list(map(lambda x: x.replace(" ", ""), ann_keys))
        if "medaka_version" in dict(variants.metadata):
            is_medaka = True
        else:
            is_medaka = False


        for variant in variants:
            ann_values = variant.INFO['ANN'][0].split("|")
            d = dict(zip(ann_keys, ann_values))
            if is_medaka:
                line = [
                    variant.CHROM,
                    variant.POS,
                    f"{variant.REF}:{sum(variant.INFO['SR'][:2])}",
                    f"{variant.ALT[0]}:{sum(variant.INFO['SR'][-2:])}",
                    variant.INFO['DP'],
                    d['Annotation'].replace("_", " "),
                    d['Annotation_Impact'],
                    d['Gene_Name'],
                    d['HGVS.c'],
                    d['HGVS.p']
                ]
            else:
                line = [
                    variant.CHROM,
                    variant.POS,
                    f"{variant.REF}:{variant.INFO['AC'][0]}",
                    f"{variant.ALT[0]}:{variant.INFO['AC'][1]}",
                    variant.INFO['DP'],
                    d['Annotation'].replace("_"," "),
                    d['Annotation_Impact'],
                    d['Gene_Name'],
                    d['HGVS.c'],
                    d['HGVS.p']
                ]
            line = map(lambda x: str(x), line)
            outf.write("\t".join(line) + "\n")


if __name__ == '__main__':
    if len(argv) != 3:
        raise Exception(usage)
    else:
        write_annotate_result(argv[1], argv[2])
