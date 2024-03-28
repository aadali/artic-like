import gzip
import argparse

from typing import List, Dict, Tuple
from collections import defaultdict

import pandas as pd
from Bio import SeqIO


# type indelInfo = Dict[str, Tuple[List[str], List[int]]]
# type region = Dict[str, List[int]]


def get_args():
    parser = argparse.ArgumentParser(
        prog="Get snv depends on the depth of each bases at each position, then "
        "Merge them and the indel of clair3 into one file."
    )
    # parser.add_argument('--dp', type=int, required=True, help="The min depth threshold for a snv")
    parser.add_argument(
        "--snv_freq",
        type=float,
        required=True,
        help="The min allele freq of alt base for a snv",
    )
    parser.add_argument(
        "--outfile",
        type=str,
        required=True,
        help="Output the variants into this vcf file",
    )
    # parser.add_argument("--low_cover_region", type=str, required=True, help="")
    parser.add_argument("reference", type=str, help="The reference fasta")
    parser.add_argument(
        "low_cover_bed", type=str, help="The low coverage region bed file"
    )
    parser.add_argument(
        "ppebd",
        type=str,
        help="A file contains 10 cols data which are "
        "contig,pos,ref,dp,a_dp,c_dp,g_dp,t_dp,del_dp,ins_dp,"
        "but only the first 8 cols will be used",
    )
    # parser.add_argument("pbd", type=str, help="Depth of per position, a file from samtools depth -a -J")
    parser.add_argument("indel", type=str, help="The indel vcf file")
    args = parser.parse_args()
    return args


def get_indel_info(indel: str) -> Dict[str, Tuple[List[str], List[int]]]:
    with (
        gzip.open(indel, "rt") if indel.endswith(".gz") else open(indel, "r")
    ) as infile:
        d = {}
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.strip("\n").split("\t")
            contig, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            indel_region = [x for x in range(int(pos), int(pos) + len(ref))]
            d[f"{contig}@@{pos}"] = (fields, indel_region)
    return d


def dont_check_snv_region(
    low_cover_bed: str, indel_info: Dict[str, Tuple[List[str], List[int]]]
) -> Dict[str, List[int]]:
    """
    merge low cover region and positions where indel happened
    :param low_cover_bed: low cover bed file
    :param indel_info: a dict, {contig}@@{position} as key,
    (fields list split by indel variant line in vcf, positions affected by indel) as value
    :return: Region
    """
    no_snv_region = defaultdict(list)
    with open(low_cover_bed, "r") as infile:
        for line in infile:
            contig, start, end = line.strip("\n").split("\t")[:3]
            start, end = int(start), int(end)
            no_snv_region[contig] += list(range(start + 1, end + 1))

    for contig_pos in indel_info:
        contig = contig_pos.split("@@")[0]
        no_snv_region[contig] += indel_info[contig_pos][1]
    return no_snv_region


def merge_snv_and_indel(
    reference: str,
    ppebd: str,
    low_cover_bed: str,
    indel_info: Dict[str, Tuple[List[str], List[int]]],
    vcf_header: str,
    outfile: str,
    min_snv_freq: float,
):
    """
    1. get region of indel and low coverage
    2. detecting snv from ppebd(per position each base dp) depends on the depth of each base at per position,
       but ignore the position in region of indel position and low coverage
    3. using min_snv_freq to filter the snv and merge the snvs and indel.
    :param reference: the reference file path
    :param ppebd: a file, depth of each base at per postion
    :param low_cover_bed: a bed file, low coverage
    :param indel_info: type indelInfo, only contains indel variants in this file, generated from calir3 and self-script
    :param vcf_header: the header of indel vcf info
    :param outfile: output file
    :param min_snv_freq: min snv freq
    :return: None
    """
    ppebd_df = pd.read_csv(ppebd, sep="\t", header=None)
    ppebd_df.columns = ["contig", "pos", "a", "c", "g", "t", "n", "del", "ins"]

    ref_records = SeqIO.parse(open(reference, "r"), format="fasta")
    ref_dict = {"contig": [], "pos": [], "ref": []}
    for record in ref_records:
        ref_seq = list(str(record.seq).upper())
        contig = [record.id] * len(ref_seq)
        pos = list(range(1, len(ref_seq) + 1))
        ref_dict["contig"] += contig
        ref_dict["pos"] += pos
        ref_dict["ref"] += ref_seq
    ref_df = pd.DataFrame.from_dict(ref_dict)
    # there is no ref_base info in the ppebd file outputed by igvtools,
    # so we need reference fasta and add ref_base column into the ppebd
    # ref_base will be used to check snv later
    ppebd_df = pd.merge(ref_df, ppebd_df, how="inner", on=["contig", "pos"])

    no_snv_region = dont_check_snv_region(low_cover_bed, indel_info)
    # print(no_snv_region)

    variants = [indel[0] for indel in list(indel_info.values())]

    for contig, ctg_no_snv_region in no_snv_region.items():
        sub_df = ppebd_df.query(
            "(contig == @contig) & (pos not in @ctg_no_snv_region)"
        ).reset_index(drop=True)
        freq_df = sub_df[list("acgt")]
        criteria1 = (
            freq_df.T.idxmax() != sub_df["ref"].str.lower()
        )  # if the most dp of base is not ref_base, it may be a mutation
        candidate_snvs: pd.DataFrame = sub_df.loc[criteria1, :].reset_index(drop=True)
        acgt = candidate_snvs[list("acgt")]
        candidate_snvs.insert(len(candidate_snvs.columns), "max_base", acgt.T.idxmax())
        candidate_snvs.insert(
            len(candidate_snvs.columns),
            "dp",
            candidate_snvs["a"]
            + candidate_snvs["t"]
            + candidate_snvs["g"]
            + candidate_snvs["c"],
        )
        for idx, var in candidate_snvs.iterrows():
            if var.dp <= 0:
                continue
            ref_dp = int(var[var.ref.lower()])
            alt_dp = int(var[var.max_base])
            """
            Calculate the percentage of alt_base in sum(ref_base, alt_base, the other two bases dp) and treat it as AF.
            Cause DP is equal to (A_DP + C_DP + G_DP + T_DP + DEL_DP) in general and the inaccuracy of ont data will 
            generate more noise such as other bases besides ref base and alt base. This make AF less than it's actual 
            value.
            Only think about binary mutation!
            """
            af = alt_dp / var.dp

            # sometimes, there are some N base in reference
            if var.ref in ["n", "N"] or af < min_snv_freq:
                continue
            ad = f"{ref_dp},{alt_dp}"
            variants.append(
                [
                    var.contig,
                    str(var.pos),
                    ".",
                    var.ref,
                    var.max_base.upper(),
                    "500",
                    "PASS",
                    "Script",
                    "DP:AD:AF",
                    f"{int(var.dp)}:{ad}:{af}",
                ]
            )
    variants.sort(key=lambda x: (x[0], int(x[1])))
    variants = map(lambda x: "\t".join(x), variants)
    with open(outfile, "w", encoding="utf-8") as outfile:
        contents = [vcf_header] + list(variants)
        outfile.write("\n".join(contents) + "\n")


def main():
    args = get_args()
    snv_freq = args.snv_freq
    outfile = args.outfile
    ppebd = args.ppebd
    reference = args.reference
    low_cover_bed = args.low_cover_bed
    indel = args.indel
    indel_info = get_indel_info(indel)
    vcf_header = ""
    with (
        gzip.open(indel, "rt") if indel.endswith(".gz") else open(indel, "r")
    ) as infile:
        for line in infile:
            if not line.startswith("#"):
                break
            vcf_header += line

    merge_snv_and_indel(
        reference=reference,
        ppebd=ppebd,
        low_cover_bed=low_cover_bed,
        indel_info=indel_info,
        vcf_header=vcf_header.strip("\n"),
        outfile=outfile,
        min_snv_freq=snv_freq,
    )


if __name__ == "__main__":
    main()
