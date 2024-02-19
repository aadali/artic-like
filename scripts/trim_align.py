import sys

from sys import argv
from os import path
from copy import copy
from collections import defaultdict

import pysam
import pandas as pd


"""
jc1200 primers infomation
if the primer_bed in config.yaml == "jc1200", the jc1200 primer will be used. Otherwise, users should specify the primer 
or not.

primer_bed info will be used to trim the primer of the amplicon
"""


def get_jc1200_bed():
    return pd.DataFrame.from_dict(
        {
            'ref': {0: 'NC_045512.2', 1: 'NC_045512.2', 2: 'NC_045512.2', 3: 'NC_045512.2', 4: 'NC_045512.2',
                    5: 'NC_045512.2',
                    6: 'NC_045512.2', 7: 'NC_045512.2', 8: 'NC_045512.2', 9: 'NC_045512.2', 10: 'NC_045512.2',
                    11: 'NC_045512.2', 12: 'NC_045512.2', 13: 'NC_045512.2', 14: 'NC_045512.2', 15: 'NC_045512.2',
                    16: 'NC_045512.2', 17: 'NC_045512.2', 18: 'NC_045512.2', 19: 'NC_045512.2', 20: 'NC_045512.2',
                    21: 'NC_045512.2', 22: 'NC_045512.2', 23: 'NC_045512.2', 24: 'NC_045512.2', 25: 'NC_045512.2',
                    26: 'NC_045512.2', 27: 'NC_045512.2', 28: 'NC_045512.2', 29: 'NC_045512.2', 30: 'NC_045512.2',
                    31: 'NC_045512.2', 32: 'NC_045512.2', 33: 'NC_045512.2', 34: 'NC_045512.2', 35: 'NC_045512.2',
                    36: 'NC_045512.2', 37: 'NC_045512.2', 38: 'NC_045512.2', 39: 'NC_045512.2', 40: 'NC_045512.2',
                    41: 'NC_045512.2', 42: 'NC_045512.2', 43: 'NC_045512.2', 44: 'NC_045512.2', 45: 'NC_045512.2',
                    46: 'NC_045512.2', 47: 'NC_045512.2', 48: 'NC_045512.2', 49: 'NC_045512.2', 50: 'NC_045512.2',
                    51: 'NC_045512.2', 52: 'NC_045512.2', 53: 'NC_045512.2', 54: 'NC_045512.2', 55: 'NC_045512.2',
                    56: 'NC_045512.2', 57: 'NC_045512.2'},
            'start': {0: 30, 1: 1183, 2: 1100, 3: 2244, 4: 2153, 5: 3235, 6: 3144, 7: 4240, 8: 4167, 9: 5337, 10: 5257,
                      11: 6358,
                      12: 6283, 13: 7379, 14: 7298, 15: 8363, 16: 8253, 17: 9378, 18: 9303, 19: 10429, 20: 10343,
                      21: 11447,
                      22: 11372, 23: 12538, 24: 12450, 25: 13599, 26: 13509, 27: 14619, 28: 14540, 29: 15713, 30: 15608,
                      31: 16698, 32: 16624, 33: 17732, 34: 17622, 35: 18684, 36: 18596, 37: 19655, 38: 19574, 39: 20676,
                      40: 20553, 41: 21613, 42: 21532, 43: 22590, 44: 22511, 45: 23609, 46: 23518, 47: 24714, 48: 24633,
                      49: 25768, 50: 25690, 51: 26835, 52: 26744, 53: 27872, 54: 27784, 55: 28985, 56: 28677,
                      57: 29768},
            'end': {0: 54, 1: 1205, 2: 1128, 3: 2266, 4: 2179, 5: 3257, 6: 3166, 7: 4262, 8: 4189, 9: 5359, 10: 5286,
                    11: 6380,
                    12: 6307, 13: 7401, 14: 7328, 15: 8385, 16: 8282, 17: 9400, 18: 9327, 19: 10451, 20: 10370,
                    21: 11469,
                    22: 11394, 23: 12560, 24: 12473, 25: 13621, 26: 13532, 27: 14641, 28: 14568, 29: 15735, 30: 15634,
                    31: 16720, 32: 16647, 33: 17754, 34: 17649, 35: 18706, 36: 18618, 37: 19678, 38: 19604, 39: 20698,
                    40: 20581, 41: 21648, 42: 21562, 43: 22612, 44: 22537, 45: 23631, 46: 23544, 47: 24736, 48: 24658,
                    49: 25790, 50: 25712, 51: 26857, 52: 26766, 53: 27894, 54: 27808, 55: 29007, 56: 28699, 57: 29790},
            'name': {0: 'SARSCoV-1200_1_LEFT', 1: 'SARSCoV-1200_1_RIGHT', 2: 'SARSCoV-1200_2_LEFT',
                     3: 'SARSCoV-1200_2_RIGHT',
                     4: 'SARSCoV-1200_3_LEFT', 5: 'SARSCoV-1200_3_RIGHT', 6: 'SARSCoV-1200_4_LEFT',
                     7: 'SARSCoV-1200_4_RIGHT',
                     8: 'SARSCoV-1200_5_LEFT', 9: 'SARSCoV-1200_5_RIGHT', 10: 'SARSCoV-1200_6_LEFT',
                     11: 'SARSCoV-1200_6_RIGHT',
                     12: 'SARSCoV-1200_7_LEFT', 13: 'SARSCoV-1200_7_RIGHT', 14: 'SARSCoV-1200_8_LEFT',
                     15: 'SARSCoV-1200_8_RIGHT', 16: 'SARSCoV-1200_9_LEFT', 17: 'SARSCoV-1200_9_RIGHT',
                     18: 'SARSCoV-1200_10_LEFT', 19: 'SARSCoV-1200_10_RIGHT', 20: 'SARSCoV-1200_11_LEFT',
                     21: 'SARSCoV-1200_11_RIGHT', 22: 'SARSCoV-1200_12_LEFT', 23: 'SARSCoV-1200_12_RIGHT',
                     24: 'SARSCoV-1200_13_LEFT', 25: 'SARSCoV-1200_13_RIGHT', 26: 'SARSCoV-1200_14_LEFT',
                     27: 'SARSCoV-1200_14_RIGHT', 28: 'SARSCoV-1200_15_LEFT', 29: 'SARSCoV-1200_15_RIGHT',
                     30: 'SARSCoV-1200_16_LEFT', 31: 'SARSCoV-1200_16_RIGHT', 32: 'SARSCoV-1200_17_LEFT',
                     33: 'SARSCoV-1200_17_RIGHT', 34: 'SARSCoV-1200_18_LEFT', 35: 'SARSCoV-1200_18_RIGHT',
                     36: 'SARSCoV-1200_19_LEFT', 37: 'SARSCoV-1200_19_RIGHT', 38: 'SARSCoV-1200_20_LEFT',
                     39: 'SARSCoV-1200_20_RIGHT', 40: 'SARSCoV-1200_21_LEFT', 41: 'SARSCoV-1200_21_RIGHT',
                     42: 'SARSCoV-1200_22_LEFT', 43: 'SARSCoV-1200_22_RIGHT', 44: 'SARSCoV-1200_23_LEFT',
                     45: 'SARSCoV-1200_23_RIGHT', 46: 'SARSCoV-1200_24_LEFT', 47: 'SARSCoV-1200_24_RIGHT',
                     48: 'SARSCoV-1200_25_LEFT', 49: 'SARSCoV-1200_25_RIGHT', 50: 'SARSCoV-1200_26_LEFT',
                     51: 'SARSCoV-1200_26_RIGHT', 52: 'SARSCoV-1200_27_LEFT', 53: 'SARSCoV-1200_27_RIGHT',
                     54: 'SARSCoV-1200_28_LEFT', 55: 'SARSCoV-1200_28_RIGHT', 56: 'SARSCoV-1200_29_LEFT',
                     57: 'SARSCoV-1200_29_RIGHT'}}

    )


consumes_ref = [1, 0, 1, 1, 0, 0, 0, 1, 1]
consumes_query = [1, 1, 0, 0, 1, 0, 0, 1, 1]

"""
cigars = {"M": 0, "I": 1, "D": 2,
          "N": 3, "S": 4, "H": 5,
          "P": 6, "=": 7, "X": 8, "B": 9}
"""


class Primer(object):
    def __init__(self, ref, start, end, name):
        self.ref = ref
        self.start = start
        self.end = end
        self.name = name

    def __repr__(self):
        return self.name


class PrimersAmplicon(object):
    def __init__(self, insert_start, insert_end, amp_start, amp_end):
        self.insert_start = insert_start
        self.insert_end = insert_end
        self.amp_start = amp_start
        self.amp_end = amp_end
        self.insert_len = insert_end - insert_start

    def __repr__(self):
        return f"{self.insert_start}-{self.insert_end}"


class Primers(object):
    def __init__(self, primer_bed):
        self.primer_bed = primer_bed
        self.primers = self.get_primers()
        # self.check_primers()

    def get_primers(self):
        if self.primer_bed == "jc1200":
            sys.stdout.write("\n\nyour choice is JC 1200 primers\n\n")
            primer_df = get_jc1200_bed()
        else:
            primer_df = pd.read_csv(self.primer_bed,
                                    sep="\t",
                                    header=None,
                                    usecols=range(4),
                                    names=['ref', 'start', 'end', 'name'],
                                    dtype={'ref': str, 'start': int, 'end': int, 'name': str})
        self.most_left = min(primer_df['start'])
        self.most_right = max(primer_df['end'])

        primers = defaultdict(dict)
        for primer in primer_df.itertuples(index=False):
            ref, start, end, name = primer
            if "_LEFT" in name:
                paired_primers_name = name.split("_LEFT")[0]
                if "left" not in primers[paired_primers_name]:
                    left_primer = Primer(ref=ref, start=start, end=end, name=name)
                    primers[paired_primers_name]['left'] = left_primer
                else:
                    raise Exception(f"Two or more left primer named {paired_primers_name} found")
            elif "_RIGHT" in name:
                paired_primers_name = name.split("_RIGHT")[0]
                if "right" not in primers[paired_primers_name]:
                    right_primer = Primer(ref=ref, start=start, end=end, name=name)
                    primers[paired_primers_name]['right'] = right_primer
                else:
                    raise Exception(f"Two or more right primer named {paired_primers_name} found")
            else:
                raise Exception(f"Couldn't detect \"_LEFT\" or \"_RIGHT\" in primer name\n")

        for paired_primers_name, paired_primers in primers.items():
            amplicon = PrimersAmplicon(
                insert_start=paired_primers['left'].end,
                insert_end=paired_primers['right'].start,
                amp_start=paired_primers['left'].start,
                amp_end=paired_primers['right'].end
            )
            primers[paired_primers_name]['amplicon'] = amplicon

        return primers

    # def check_primers(self):
    #     for paired_primers_name, paired_primers in self.primers.items():
    #         if len(paired_primers) != 2:
    #             raise Exception(f"More or fewer primers found: {paired_primers_name}")
    #         left_primer, right_primer = paired_primers['left'], paired_primers['right']
    #         if (left_primer['name'].split("_LEFT")[0] != paired_primers_name or
    #                 right_primer['name'].split("_RIGHT")[0] != paired_primers_name):
    #             raise Exception(f"Not paired primers found:{paired_primers_name}")


def find_closest_amplicon(segment, primers):
    assert isinstance(primers, Primers)
    assert isinstance(segment, pysam.AlignedSegment)
    ref_start, ref_end = segment.reference_start, segment.reference_end
    candidate_primers = []
    for paired_primers_name, paired_primers in primers.primers.items():
        amplicon = paired_primers['amplicon']
        if ref_start <= amplicon.insert_start <= ref_end <= amplicon.insert_end:
            """
            situation 1:
                        insert_start-------------------------------insert_end
            ref_start---------------------------------ref_end
            """
            candidate_primers.append([paired_primers, abs(ref_end - amplicon.insert_start)])

        elif ref_start <= amplicon.insert_start <= amplicon.insert_end <= ref_end:
            """
            situation 2:
                        insert_start-------------------------------insert_end
            ref_start--------------------------------------------------------------ref_end
            """
            candidate_primers.append([paired_primers, abs(amplicon.insert_end - amplicon.insert_start)])

        elif amplicon.insert_start <= ref_start <= amplicon.insert_end <= ref_end:
            """
            situation 3:
            insert_start------------------------------insert_end
                        ref_start----------------------------------ref_end
            """
            candidate_primers.append([paired_primers, abs(amplicon.insert_end - ref_start)])

        elif amplicon.insert_start <= ref_start <= ref_end <= amplicon.insert_end:
            """
            situation 4:
            insert_start--------------------------------------------------------------insert_end
                        ref_start-------------------------------ref_end

            """
            candidate_primers.append([paired_primers, abs(ref_start - ref_end)])

        else:
            continue

    if candidate_primers:
        return sorted(candidate_primers, key=lambda x: x[1], reverse=True)[0][0]
    else:
        return


def trim_left(segment, primers=None, trim_bases=None, most_left=False, supplementary=False):
    assert isinstance(segment, pysam.AlignedSegment)
    segment = copy(segment)
    left_soft = 0
    offset = 0
    left_cigar_modify = []
    ref_start, ref_end = segment.reference_start, segment.reference_end
    cigars = segment.cigartuples
    first_cigar = cigars[0]
    if primers is not None:
        #  for the first amplicon's left primer, trim from the primer start
        left_end = primers['left'].start if most_left else primers['left'].end
    elif trim_bases is not None:
        assert isinstance(trim_bases, int)
        left_end = ref_start + trim_bases

    else:
        raise Exception("only one params shoud be spicified in {\"primers\", \"trim_bases\"}")

    if ref_start < left_end:
        segment.reference_start = left_end
        total_offset = left_end - ref_start
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
                        segment.reference_start += offset - total_offset
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
                segment.query_sequence = segment.query_sequence[left_soft:]
            else:
                cigars.append((4, left_soft))
            cigars.reverse()
            segment.cigartuples = cigars
            return segment
        cigars.append((4, left_soft))
        cigars.reverse()
    segment.cigartuples = cigars
    return segment


def trim_right(segment, primers=None, trim_bases=None, most_right=False, supplementary=False):
    assert isinstance(segment, pysam.AlignedSegment)
    segment = copy(segment)
    right_cigar_modify = []
    right_soft = 0
    offset = 0
    cigars = segment.cigartuples
    last_cigar = cigars[-1]
    ref_start, ref_end = segment.reference_start, segment.reference_end
    if primers is not None:
        #  for the last amplicon's right primer, trim from the primer end
        right_start = primers['right'].end if most_right else primers['right'].start
    elif trim_bases is not None:
        assert isinstance(trim_bases, int)
        right_start = ref_end - trim_bases
    else:
        raise Exception("only one params shoud be spicified in {\"primers\", \"trim_bases\"}")

    if ref_end > right_start:
        total_offset = ref_end - right_start
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
            query_sequence = segment.query_sequence
            if last_cigar[0] == 5:
                # for supplementary align, if the last cigar is "H",
                # so right_soft + the original  "H" number bases need to be hard clipped,
                # sequence filed also need update
                segment.query_sequence = query_sequence[:-right_soft] if right_soft > 0 else query_sequence
                cigars.append((5, right_soft + last_cigar[1]))
            else:
                cigars.append((4, right_soft))
            segment.cigartuples = cigars
            return segment
        cigars.append((4, right_soft))
    segment.cigartuples = cigars
    return segment


def trim_bam(inbam_fp, outbam_fp, primer_bed, trim_log, min_overlap, trim_bases):
    inbam = pysam.AlignmentFile(inbam_fp, "rb")
    outbam = pysam.AlignmentFile(outbam_fp, "wb", template=inbam)
    trim_log_fh = open(trim_log, 'w', encoding='utf-8')
    primer_base_name = path.basename(primer_bed)
    no_primer_header = "\t".join(['ReadName', "RefStart", "RefEnd", "AlignLen",
                                  "RefStartTrimmed", "RefEndTrimmed", "AlignLenTrimmed"])
    header = "\t".join([
        "ReadName",
        "RefStart",
        "RefEnd",
        "AlignLen",
        "PairedPrimerName",
        "LeftPrimer",
        # "LeftPrimerStart",
        "LeftPrimerEnd",
        "RightPrimer",
        "RightPrimerStart",
        # "RightPrimerEnd",
        "RefStartTrimmed",
        "RefEndTrimmed",
        "AlignLenTrimmed",
        "InsertLen"
    ])
    if primer_base_name == ".noprimer.txt":
        trim_bases = int(trim_bases)
        trim_log_fh.write(no_primer_header + '\n')
        for segment in inbam:
            if segment.is_unmapped or segment.is_supplementary:
                continue
            new_segment = trim_left(segment, trim_bases=trim_bases)
            trimmed_segment = trim_right(new_segment, trim_bases=trim_bases)
            if trimmed_segment.reference_length < int(min_overlap):
                continue
            trim_log_fh.write(f"{segment.qname}\t"
                              f"{segment.reference_start}\t"
                              f"{segment.reference_end}\t"
                              f"{segment.reference_length}\t"
                              f"{trimmed_segment.reference_start}\t"
                              f"{trimmed_segment.reference_end}\t"
                              f"{trimmed_segment.reference_length}\n")
            outbam.write(trimmed_segment)

    else:
        primers = Primers(primer_bed)
        trim_log_fh.write(header + '\n')
        for segment in inbam:
            if segment.is_unmapped or segment.is_supplementary:
                continue
            paired_primers = find_closest_amplicon(segment, primers=primers)
            if paired_primers is None:
                continue
            else:
                most_left = paired_primers['left'].start == primers.most_left
                mosrt_right = paired_primers['right'].end == primers.most_right
                new_segment = trim_left(segment, primers=paired_primers, most_left=most_left)
                trimmed_segment = trim_right(new_segment, primers=paired_primers, most_right=mosrt_right)
                record = "\t".join([
                    segment.qname,
                    str(segment.reference_start),
                    str(segment.reference_end),
                    str(segment.reference_length),
                    paired_primers['left'].name.split("_LEFT")[0],
                    paired_primers['left'].name,
                    # str(paired_primers['left'].start),
                    str(paired_primers['left'].end),
                    paired_primers['right'].name,
                    str(paired_primers['right'].start),
                    # str(paired_primers['right'].end),
                    str(trimmed_segment.reference_start),
                    str(trimmed_segment.reference_end),
                    str(trimmed_segment.reference_length),
                    str(paired_primers['amplicon'].insert_len)
                ]) + "\n"
                if trimmed_segment.reference_length < int(min_overlap):
                    continue
                trim_log_fh.write(record)
                outbam.write(trimmed_segment)


if __name__ == '__main__':
    usage = f"usage: {path.basename(__file__)} <input_bam> <output_bam> <primer_bed> <trim_log> <min_overlap> <trim_bases>"
    paras_num = len(argv)
    if paras_num != 7:
        raise Exception(usage)
    else:
        trim_bam(*argv[1:])
