from sys import argv
from copy import copy
import pandas as pd
from collections import defaultdict
from os import path
import pysam

consumes_ref = [1, 0, 1, 1, 0, 0, 0, 1, 1]
consumes_query = [1, 1, 0, 0, 1, 0, 0, 1, 1]
usage = f"usage: {path.basename(__file__)} <input_bam> <output_bam> <primer_bed> <trim_log>"


class Primer(object):
    def __init__(self, ref, start, end, name):
        self.ref = ref
        self.start = start
        self.end = end
        self.name = name


class PrimersAmplicon(object):
    def __init__(self, insert_start, insert_end, amp_start, amp_end):
        self.insert_start = insert_start
        self.insert_end = insert_end
        self.amp_start = amp_start
        self.amp_end = amp_end
        self.insert_len = insert_end - insert_start


class Primers(object):
    def __init__(self, primer_bed):
        self.primer_bed = primer_bed
        self.primers = self.get_primers()
        # self.check_primers()

    def get_primers(self):
        primer_df = pd.read_csv(self.primer_bed,
                                sep="\t",
                                header=None,
                                usecols=range(4),
                                names=['ref', 'start', 'end', 'name'],
                                dtype={'ref': str, 'start': int, 'end': int, 'name': str})

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
                raise Exception(f"Couldn't detect \"_LEFT\" or \"_RIGHT\" in primer name\n{usage}")

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
        if amplicon.insert_start >= ref_start and amplicon.insert_end >= ref_end:
            """
            situation 1:
                        insert_start-------------------------------insert_end
            ref_start---------------------------------ref_end
            """
            candidate_primers.append([paired_primers, ref_end - amplicon.insert_start])

        elif amplicon.insert_start >= ref_start and amplicon.insert_end <= ref_end:
            """
            situation 2:
                        insert_start-------------------------------insert_end
            ref_start--------------------------------------------------------------ref_end
            """
            candidate_primers.append([paired_primers, amplicon.insert_end - amplicon.insert_start])

        elif amplicon.insert_start <= ref_start and amplicon.insert_end <= ref_end:
            """
            situation 3:
            insert_start------------------------------insert_end
                        ref_start----------------------------------ref_end
            """
            candidate_primers.append([paired_primers, amplicon.insert_end - ref_start])

        elif amplicon.insert_start <= ref_start and amplicon.insert_end >= ref_end:
            """
            situation 4:
            insert_start--------------------------------------------------------------insert_end
                        ref_start-------------------------------ref_end
            
            """
            candidate_primers.append([paired_primers, ref_start - ref_end])

        else:
            continue

    if candidate_primers:
        return max(candidate_primers, key=lambda x: x[1])[0]
    else:
        return


def trim_left(segment, primers, supplementary=False):
    assert isinstance(segment, pysam.AlignedSegment)
    segment = copy(segment)
    left_soft = 0
    offset = 0
    left_cigar_modify = []
    ref_start, ref_end = segment.reference_start, segment.reference_end
    cigars = segment.cigartuples
    first_cigar = cigars[0]
    left_primer = primers['left']

    if ref_start < left_primer.end:
        segment.reference_start = left_primer.end
        total_offset = left_primer.end - ref_start
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


def trim_right(segment, primers, supplementary=False):
    assert isinstance(segment, pysam.AlignedSegment)
    segment = copy(segment)
    right_cigar_modify = []
    right_soft = 0
    offset = 0
    cigars = segment.cigartuples
    last_cigar = cigars[-1]
    ref_start, ref_end = segment.reference_start, segment.reference_end
    right_primer = primers['right']

    if ref_end > right_primer.start:
        total_offset = ref_end - right_primer.start
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


def trim_bam(inbam_fp, outbam_fp, primer_bed, trim_log):
    inbam = pysam.AlignmentFile(inbam_fp, "rb")
    outbam = pysam.AlignmentFile(outbam_fp, "wb", template=inbam)
    trim_log_fh = open(trim_log, 'w', encoding='utf-8')
    primer_base_name = path.basename(primer_bed)
    no_primer_header = "\t".join(['ReadName', "RefStart", "RefEnd", "AlignLen"])
    header = "\t".join([
        "ReadName",
        "RefStart",
        "RefEnd",
        "AlignLen",
        "PairedPrimerName",
        "LeftPrimer",
        "LeftPrimerStart",
        "LeftPrimerEnd",
        "RightPrimer",
        "RightPrimerStart",
        "RightPrimerEnd",
        "RefStartTrimmed",
        "RefEndTrimmed",
        "AlignLenTrimmed",
        "InsertLen"
    ])
    if primer_base_name == ".noprimer.txt":
        trim_log_fh.write(no_primer_header + '\n')
        for segment in inbam:
            if segment.is_unmapped or segment.is_supplementary:
                continue
            trim_log_fh.write(f"{segment.qname}\t"
                              f"{str(segment.reference_start)}\t"
                              f"{str(segment.reference_end)}\t"
                              f"{segment.reference_length}\n")
            outbam.write(segment)

    else:
        primers = Primers(primer_bed)
        trim_log_fh.write(header + '\n')
        for segment in inbam:
            if segment.is_unmapped or segment.is_supplementary:
                continue
            paired_primers = find_closest_amplicon(segment, primers)
            if paired_primers is None:
                continue
            else:
                new_segment = trim_left(segment, paired_primers)
                trimmed_segment = trim_right(new_segment, paired_primers)
                record = "\t".join([
                    segment.qname,
                    str(segment.reference_start),
                    str(segment.reference_end),
                    str(segment.reference_length),
                    paired_primers['left'].name.split("_LEFT")[0],
                    paired_primers['left'].name,
                    str(paired_primers['left'].start),
                    str(paired_primers['left'].end),
                    paired_primers['right'].name,
                    str(paired_primers['right'].start),
                    str(paired_primers['right'].end),
                    str(trimmed_segment.reference_start),
                    str(trimmed_segment.reference_end),
                    str(trimmed_segment.reference_length),
                    str(paired_primers['amplicon'].insert_len)
                ]) + "\n"
                if trimmed_segment.reference_length < 100:
                    continue
                trim_log_fh.write(record)
                outbam.write(trimmed_segment)


if __name__ == '__main__':
    paras_num = len(argv)
    if paras_num != 5:
        raise Exception(usage)
    else:
        trim_bam(*argv[1:])
