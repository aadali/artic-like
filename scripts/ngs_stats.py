import re
from sys import argv

fastp_log = argv[1]
bam_stats = argv[2]
outfile = argv[3]

with open(fastp_log) as infile1:
    filt_result = infile1.read()
    reads_info = re.findall("total reads: (\d+)", filt_result)
    bases_info = re.findall("total bases: (\d+)", filt_result)
    # q20_info = re.findall("Q20 bases: ([\d\.\(\)%]+)", filt_result)
    q30_info = re.findall("Q30 bases: ([\d\.\(\)%]+)", filt_result)

with open(bam_stats) as infile2:
    bam_stats_result = infile2.read()
    mapped_reads = re.findall("reads mapped:\s+(\d+)", bam_stats_result)
    mapped_bases = re.findall("bases mapped:\s+(\d+)", bam_stats_result)
    reads_paired_mapped = re.findall("reads mapped and paired:\s+(\d+)", bam_stats_result)
    reads_paired_mapped = reads_paired_mapped if reads_paired_mapped else ['0']

with open(outfile, 'w') as outfile:
    if len(reads_info) == 4:
        outfile.write(
            f"raw reads1 reads\t{reads_info[0]}\n"
            f"clean reads1 reads\t{reads_info[2]}\n"
            f"raw reads1 bases\t{bases_info[0]}\n"
            f"clean reads1 bases\t{bases_info[2]}\n"
            f"raw reads2 reads\t{reads_info[1]}\n"
            f"clean reads2 reads\t{reads_info[3]}\n"
            f"raw reads2 bases\t{bases_info[1]}\n"
            f"clean reads2 bases\t{bases_info[3]}\n"
            f"raw reads1 Q30\t{q30_info[0]}\n"
            f"clean reads1 Q30\t{q30_info[2]}\n"
            f"raw reads2 Q30\t{q30_info[1]}\n"
            f"clean reads2 Q30\t{q30_info[3]}\n"
            f"mapped_reads\t{mapped_reads[0]}\n"
            f"paired_mapped_reads\t{reads_paired_mapped[0]}\n"
            f"mapped_bases\t{mapped_bases[0]}\n"
        )
    else:
        outfile.write(
            f"raw reads1 reads\t{reads_info[0]}\n"
            f"clean reads1 reads\t{reads_info[1]}\n"
            f"raw reads1 bases\t{bases_info[0]}\n"
            f"clean reads1 bases\t{bases_info[1]}\n"
            f"raw reads1 Q30\t{q30_info[0]}\n"
            f"clean reads1 Q30\t{q30_info[1]}\n"
            f"mapped_reads\t{mapped_reads[0]}\n"
            f"paired_mapped_reads\t{reads_paired_mapped[0]}\n"
            f"mapped_bases\t{mapped_bases[0]}\n"
        )

