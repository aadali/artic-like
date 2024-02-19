import argparse
from argparse import Namespace

"""
从原始的freebayes输出vcf文件中过滤掉如下record；
    1. 有多个ALT碱基的record
    2. REF深度与ALT深度之和小于指定参数 min_dp的record
    3. REF深度与ALT深度之和中ALT深度小于指定参数 min_snv_freq的record
注：计算深度时不是使用软件输出的DP，而是使用REF深度与ALT深度之和，对于不是二元突变的位置（所以要去掉非二元突变）
    或者低质量碱基太多（freebayes命令中过滤掉低于指定质量值的碱基）的位置上来说，DP和REF、ALT深度之和会有一定差距。
"""


def get_args() -> Namespace:
    parser = argparse.ArgumentParser("Filter the output vcf of freebayes depends on some custom criteria")
    parser.add_argument("--infile", type=str, required=True,
                        help="the input vcf file")
    parser.add_argument("--outfile", type=str, required=True,
                        help="the output vcf file")
    parser.add_argument("--min_dp", type=int, default=10,
                        help="Require at least this coverage to process a site. default: 10")
    parser.add_argument("--min_alt_freq", type=float, default=0.05,
                        help="Require at least this count of observations supporting an alternate allele within a single "
                             "individual in order to evaluate the position.  default: 0.05")
    args = parser.parse_args()
    return args


def filter_record(record: str, min_dp: int, min_alt_freq: float) -> bool:
    if record.startswith("#"):
        return True
    fields = record.strip().split("\t")
    alt, fmt, sample = fields[4], fields[-2], fields[-1]
    if "," in alt:
        return False
    ad_idx = fmt.split(":").index('AD')
    ref_dp, alt_dp = list(map(lambda x: int(x), sample.split(":")[ad_idx].split(",")))
    if ref_dp + alt_dp < min_dp or alt_dp / (ref_dp + alt_dp) < min_alt_freq:
        return False
    return True


def main() -> None:
    args = get_args()
    with open(args.infile) as infile:
        contents = [line for line in infile if filter_record(line, args.min_dp, args.min_alt_freq)]
    with open(args.outfile, 'w') as outfile:
        outfile.write("".join(contents))


if __name__ == '__main__':
    main()
