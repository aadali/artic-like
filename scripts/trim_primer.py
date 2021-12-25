from os import path
from sys import argv
from collections import defaultdict

import pysam
import pandas as pd

usage = f"usage: {path.basename(__file__)} <input_bam> <output_bam> [primer_bed]"

"""
cigars = {"M": 0, "I": 1, "D": 2,
          "N": 3, "S": 4, "H": 5,
          "P": 6, "=": 7, "X": 8, "B": 9}
"""

consumes_ref = [1, 0, 1, 1, 0, 0, 0, 1, 1]
consumes_query = [1, 1, 0, 0, 1, 0, 0, 1, 1]


class Primers(object):
    def __init__(self, primer_bed):
        self.primer_bed = primer_bed
        self.primers = self.get_primers()
        self.check_primers()
        self.left_primers = [self.primers[paired_primers_name]['left']
                             for paired_primers_name in self.primers]
        self.right_primers = [self.primers[paired_primers_name]['right']
                              for paired_primers_name in self.primers]

    def get_primers(self):
        primer_df = pd.read_csv(self.primer_bed,
                                sep="\t",
                                usecols=range(4),
                                names=['ref', 'start', 'end', 'name'],
                                dtype={'ref': str, 'start': int, 'end': int, 'name': str})
        primers = defaultdict(dict)
        for primer in primer_df.itertuples(index=False):
            ref, start, end, name = primer

            if "_LEFT" in name:
                paired_primers_name = name.split("_LEFT")[0]
                if 'left' not in primers[paired_primers_name]:
                    primers[paired_primers_name]['left'] = dict(ref=ref, start=start, end=end, name=name)
                else:
                    raise Exception(f"Two or more left primer named {paired_primers_name} found")
            elif "_RIGHT" in name:
                paired_primers_name = name.split("_RIGHT")[0]
                if "right" not in primers[paired_primers_name]:
                    primers[paired_primers_name]['right'] = dict(ref=ref, start=start, end=end, name=name)
                else:
                    raise Exception(f"Two or more right primer named {paired_primers_name} found")
            else:
                raise Exception(f"Couldn't detect \"_LEFT\" or \"_RIGHT\" in primer name\n{usage}")
        return primers

    def check_primers(self):
        for paired_primer_name, paired_primers in self.primers.items():
            if len(paired_primers) != 2:
                raise Exception(f"More or fewer primers found: {paired_primer_name}")

            left_primer = paired_primers['left']
            right_primer = paired_primers['right']
            if (left_primer['name'].split("_LEFT")[0] != paired_primer_name or
                    right_primer['name'].split("_RIGHT")[0] != paired_primer_name):
                raise Exception(f"Not paired primers found: {paired_primer_name}")


def trim_left(align, left_primer):
    left_soft = 0
    offset = 0
    left_cigar_modify = []
    assert isinstance(align, pysam.AlignedSegment)
    ref_start, ref_end = align.reference_start, align.reference_end
    cigars = align.cigartuples
    if ref_start < left_primer['end']:
        align.reference_start = left_primer['end']
        total_offset = left_primer['end'] - ref_start
        for idx, (flag, length) in enumerate(cigars):
            if consumes_query[flag]:
                left_soft += length

            if consumes_ref[flag]:
                offset += length
                if offset >= total_offset:
                    if consumes_query[flag]:
                        left_soft -= (offset - total_offset)
                        left_cigar_modify = [idx, offset - total_offset]
                        # the cigar at (idx) position retain (offset-total_offset) bases

                    if flag == 2:  # flag is D, avoid situation like 30S4D, 20S4D and so on
                        left_cigar_modify = [idx, 0]
                        # the cigar at (idx) position is DEL, DEL consume ref but not query
                        align.reference_start += offset - total_offset
                    break

    if left_cigar_modify:
        idx, retain_bases = left_cigar_modify
        if retain_bases == 0:
            cigars = cigars[idx + 1:]
        else:
            cigars = cigars[idx:]
            cigars[0] = (cigars[0][0], retain_bases)
        cigars.reverse()
        cigars.append((4, left_soft))
        cigars.reverse()
    align.cigartuples = cigars
    return align


def trim_right(align, right_primer):
    assert isinstance(align, pysam.AlignedSegment)
    right_cigar_modify = []
    right_soft = 0
    offset = 0
    cigars = align.cigartuples
    ref_start, ref_end = align.reference_start, align.reference_end

    if ref_end > right_primer['start']:
        total_offset = ref_end - right_primer['start']
        for idx, (flag, length) in enumerate(cigars[::-1]):
            if consumes_query[flag]:
                right_soft += length

            if consumes_ref[flag]:
                offset += length

                if offset >= total_offset:
                    if consumes_query[flag]:
                        right_soft -= (offset - total_offset)
                        right_cigar_modify = [-idx - 1, offset - total_offset]
                    if flag == 2:  # flag is D, consuming ref
                        right_cigar_modify = [-idx - 1, 0]
                    break

    if right_cigar_modify:
        idx, retain_bases = right_cigar_modify
        cigars = cigars[:idx + 1] if idx != -1 else cigars
        if retain_bases == 0:
            cigars.pop()
        else:
            cigars[-1] = (cigars[-1][0], retain_bases)
        cigars.append((4, right_soft))
    align.cigartuples = cigars
    return align


def find_closest_primers(primers, align):
    assert isinstance(primers, Primers)
    assert isinstance(align, pysam.AlignedSegment)
    left_primer_distances = [(left_primer, abs(left_primer['start'] - align.reference_start))
                             for left_primer in primers.left_primers]
    right_primer_distances = [(right_primer, abs(right_primer['end'] - align.reference_end))
                              for right_primer in primers.right_primers]
    closest_left_primer = min(left_primer_distances, key=lambda x: x[1])[0]
    closest_right_primer = min(right_primer_distances, key=lambda x: x[1])[0]
    return closest_left_primer, closest_right_primer


def trim_bam(inbam_fp, outbam_fp, primers):
    inbam = pysam.AlignmentFile(inbam_fp, "rb")
    outbam = pysam.AlignmentFile(outbam_fp, "wb", template=inbam)
    primer_base_name = path.basename(primers)
    if primer_base_name == ".noprimer.txt":
        for align in inbam:
            if align.is_unmapped or align.is_supplementary:
                continue
            outbam.write(align)
    else:
        primers = Primers(primers)
        for align in inbam:
            left_primer, right_primer = find_closest_primers(primers, align)
            is_paired_primers = left_primer['name'].split("_LEFT")[0] == right_primer['name'].split("_RIGHT")[0]

            # drop unmapped, not paired primers and supplementary aligns
            if align.is_unmapped or not is_paired_primers or align.is_supplementary:
                continue

            left_trimmed_align = trim_left(align, left_primer)
            trimmed_align = trim_right(left_trimmed_align, right_primer)
            outbam.write(trimmed_align)


if __name__ == '__main__':
    # usage = f"usage: {path.basename(__file__)} <input_bam> <output_bam> <primer_bed>"
    paras_num = len(argv)
    if paras_num != 4:
        raise Exception(usage)
    else:
        trim_bam(argv[1], argv[2], argv[3])
