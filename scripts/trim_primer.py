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

    def get_circle_primers(self):
        """
        :return: one paired primers, the pcr production of this primer will span the head and tail of the circle genome
        or
        a empty dict
        """
        primers = self.get_primers()
        for paired_primers_name in primers:
            left_pos = primers[paired_primers_name]['left']['start']
            right_pos = primers[paired_primers_name]['right']['end']
            if right_pos < left_pos:
                return primers[paired_primers_name]
        return


def paired_supplementary(align1, align2, primers, tail_head_primers, reference_len):
    left_primer1, right_primer1 = find_closest_primers(primers, align1)
    left_primer2, right_primer2 = find_closest_primers(primers, align2)
    tail_primer, head_primer = tail_head_primers['left'], tail_head_primers['right']
    # situation 1
    situation1 = (left_primer1['start'] == tail_primer['start'] and
                  abs(align1.reference_end - reference_len) <= abs(right_primer1['end'] - align1.reference_end) and
                  right_primer2['end'] == head_primer['end'] and
                  align2.reference_start <= abs(left_primer2['start'] - align2.reference_start))

    # situation 2
    situation2 = (left_primer2['start'] == tail_primer['start'] and
                  abs(align2.reference_end - reference_len) <= abs(right_primer2['end'] - align2.reference_end) and
                  right_primer1['end'] == head_primer['end'] and
                  align1.reference_start <= abs(left_primer1['start'] - align1.reference_start))
    if situation1:
        return {'left': align1, 'right': align2}
    if situation2:
        return {'left': align2, 'right': align1}
    return None


def get_supplementary(bam_fp, primers, tail_head_primers):
    """
    get all segment that span the last base and the first base of reference.
    this segment was cloned by the tail primer (the last primer on the reference position) and
    the head primer (the first primer on the reference position)
    :param bam_fp: bam file path to be trimmed
    :param tail_head_primers: dict: {"left": tail primer, "right": head primer}
    :return: a dict containing all segments that can align to the tail reference and the head reference,
    cloned by tail primer and head primer
    """
    if tail_head_primers is None:
        return None
    bam = pysam.AlignmentFile(bam_fp, "rb")
    reference_len = bam.header.get("SQ")[0]['LN']
    d = {}
    for segment in bam:
        ref_start = segment.reference_start
        ref_end = segment.reference_end
        qname = segment.qname
        supp = segment.is_supplementary
        strand = "-" if segment.is_reverse else "+"
        if qname not in d:
            d[qname] = [dict(ref_start=ref_start,
                             ref_end=ref_end,
                             supp=supp,
                             strand=strand,
                             seg=segment)]
        else:
            d[qname].append(dict(ref_start=ref_start,
                                 ref_end=ref_end,
                                 supp=supp,
                                 strand=strand,
                                 seg=segment))

    # only reserve the segment that can only align two positions and
    # one align is primary-align,
    # the left closest is tail primer and the right closest is reference first base
    # the other align is supplementary-align,
    # the left closest is reference first base and the right closest is the head primer
    d2 = {}
    for key in d:
        if (len(d[key]) == 2) and (d[key][0]['strand'] == d[key][1]['strand']) and (
                d[key][0]['supp'] + d[key][1]['supp'] == 1):
            segment1, segment2 = d[key][0]['seg'], d[key][1]['seg']
            paired_supp_aligns = paired_supplementary(segment1, segment2, primers, tail_head_primers, reference_len)
            if paired_supp_aligns:
                d2[key] = [paired_supp_aligns['left'], paired_supp_aligns['right']]
    return d2


def trim_left(align, left_primer, supplementary=False):
    left_soft = 0
    offset = 0
    left_cigar_modify = []
    assert isinstance(align, pysam.AlignedSegment)
    ref_start, ref_end = align.reference_start, align.reference_end
    cigars = align.cigartuples
    first_cigar = cigars[0]
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
        if supplementary:
            if first_cigar[0] == 5:
                cigars.append((5, left_soft + first_cigar[1]))
                align.query_sequence = align.query_sequence[left_soft:]
            else:
                cigars.append((4, left_soft))
            cigars.reverse()
            align.cigartuples = cigars
            return align
        cigars.append((4, left_soft))
        cigars.reverse()
    align.cigartuples = cigars
    return align


def trim_right(align, right_primer, supplementary=False):
    assert isinstance(align, pysam.AlignedSegment)
    right_cigar_modify = []
    right_soft = 0
    offset = 0
    cigars = align.cigartuples
    last_cigar = cigars[-1]
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
        if supplementary:
            query_sequence = align.query_sequence
            if last_cigar[0] == 5:
                # for supplementary align, if the last cigar is "H",
                # so right_soft + the original  "H" number bases need to be hard clipped,
                # sequence filed also need update
                align.query_sequence = query_sequence[:-right_soft] if right_soft > 0 else query_sequence
                cigars.append((5, right_soft + last_cigar[1]))
            else:
                cigars.append((4, right_soft))
            align.cigartuples = cigars
            return align
        cigars.append((4, right_soft))
    align.cigartuples = cigars
    return align


def trim_supplement_align(align, candidate_supp_aligns, left_primer, right_primer):
    assert isinstance(align, pysam.AlignedSegment)
    left, right = candidate_supp_aligns[align.qname]
    if align.reference_start == left.reference_start:
        left = True
    elif align.reference_start == right.reference_start:
        left = False
    else:
        raise Exception()

    if left:  # trim left
        if not align.is_supplementary:
            new_segment = trim_left(align, left_primer)
        else:
            new_segment = trim_left(align, left_primer, supplementary=True)

    else:  # trim right
        if not align.is_supplementary:
            new_segment = trim_right(align, right_primer)
        else:
            new_segment = trim_right(align, right_primer, supplementary=True)

    return new_segment


def trim_normal_align(align, left_primer, right_primer, debug=True):
    if debug:
        return align

    new_segment = trim_left(align, left_primer)
    new_segment = trim_right(new_segment, right_primer)
    return new_segment


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
        tail_head_primers = primers.get_circle_primers()
        candidate_supp_aligns = get_supplementary(inbam_fp, primers, tail_head_primers)
        for align in inbam:
            left_primer, right_primer = find_closest_primers(primers, align)
            is_paired_primers = left_primer['name'].split("_LEFT")[0] == right_primer['name'].split("_RIGHT")[0]

            # drop unmapped, not paired primers and supplementary aligns
            if align.is_unmapped:
                continue

            if tail_head_primers is None:
                if align.is_supplementary:
                    continue
            else:
                if align.qname in candidate_supp_aligns:
                    # trim the segment depend on the align position (the tail portion or the head poftion of the reference)
                    trimed_align = trim_supplement_align(align, candidate_supp_aligns, left_primer, right_primer)
                    outbam.write(trimed_align)
                else:
                    if align.is_supplementary:
                        continue

            if is_paired_primers:
                # trim the both ends of the align
                trimed_align = trim_normal_align(align, left_primer, right_primer, debug=False)
                outbam.write(trimed_align)


if __name__ == '__main__':
    # usage = f"usage: {path.basename(__file__)} <input_bam> <output_bam> <primer_bed>"
    paras_num = len(argv)
    if paras_num != 4:
        raise Exception(usage)
    else:
        trim_bam(argv[1], argv[2], argv[3])
