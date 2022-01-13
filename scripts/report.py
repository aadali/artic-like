from os import path
import re
from datetime import datetime
from sys import argv

import pandas as pd

usage = f"{path.basename(__file__)} <SAMPLE_NAME>  <WHAT_SAMPLE> <REPORT_TEX>"
if len(argv) != 4:
    raise Exception(usage + "\n")
sample_name = argv[1]
what_sample = argv[2]
report_tex = argv[3]

# sample_name = "test001"
# what_sample = "sars-cov-2"
time = datetime.now().strftime("%Y-%m-%d %H:%M")

preamble = r"""
\documentclass[10pt]{article}
\usepackage{ctex}
\usepackage{longtable}
\usepackage{tabularx}
\usepackage{fontspec}
\usepackage{array}
\usepackage{graphicx}
\usepackage{longtable}
\usepackage{vtable}
\usepackage{geometry}
\usepackage{subfigure}
\setCJKmainfont{SourceHanSansCN-Regular}[
    Extension = .otf,
    Path = FontsPath/,
    BoldFont = SourceHanSansCN-Bold
]
\setmainfont{SourceHanSansCN-Regular}[
    Extension = .otf,
    Path = FontsPath/
]

\geometry{left=1cm, right=1cm}
\setlength{\tabcolsep}{1pt}
\renewcommand{\arraystretch}{1.3}

"""

scripts_dir = path.dirname(path.abspath(__file__))
fonts_path = path.join(scripts_dir, "../texlive/fonts")
preamble = preamble.replace("FontsPath", fonts_path)


def get_files(name):
    files = dict(
        stat_fp=path.abspath(f"{name}/stat/{name}.stat"),
        cov_fig_fp=path.abspath(f"{name}/figures/{name}.coverage.pdf"),
        annotated_var_list=path.abspath(f"{name}/variants/{name}.report.snpEff.annotate.txt"),
        unannotated_vcf_fp=path.abspath(f"{name}/variants/{name}.report.vcf"),
        consensus_fp=path.abspath(f"{name}/consensus/{name}.consensus.fasta"),
        lineage_report_fp=path.abspath(f"{name}/pangolin/{name}.lineage_report.csv"),
        nanoplot_raw_pic=path.abspath(f"{name}/nanoplot/{name}_raw.LengthvsQualityScatterPlot_dot.png"),
        nanoplot_clean_pic=path.abspath(f"{name}/nanoplot/{name}_clean.LengthvsQualityScatterPlot_dot.png"),
        nanoplot_qc_summary=path.abspath(f"{name}/nanoplot/{name}.qc.summary.txt")
    )
    return files


def escape_char(element):
    element = str(element)
    return (element.
            replace("#", "\\#").
            replace("$", "\\$").
            replace("%", "\\%").
            replace("{", "\\{").
            replace("}", "\\}").
            replace("_", "\\_").
            replace("&", "\\&")
            )


def infomations(stat_file):
    tex = []
    with open(stat_file, 'r') as infile:
        lines = infile.readlines()

    raw = re.search("(\d+)\W+(\d+)", lines[0].strip())
    raw_line_num, raw_data = raw.group(1), raw.group(2)
    raw_data = f"{int(raw_data) / 1000000:.2f}M"

    mapped = re.search("(\d+)\W+(\d+)", lines[1].strip())
    mapped_line_num, mapped_data = mapped.group(1), mapped.group(2)
    mapped_data = f"{int(mapped_data) / 1000000:.2f}M"
    tex.append("\\section{基本信息}")
    tex.append("\\begin{tabular}[htbp]{>{\\raggedright\\arraybackslash}m{5cm} >{\\raggedright\\arraybackslash}m{5cm}}")
    tex.append(f"样本名称 & {escape_char(sample_name)} \\\\")
    tex.append(f"分析时间 & {escape_char(time)} \\\\")
    tex.append(f"测序数据量 & {escape_char(raw_data)} \\\\")
    tex.append(f"测序read数 & {escape_char(raw_line_num)} \\\\")
    tex.append(f"{escape_char(what_sample)}数据量 & {escape_char(mapped_data)} \\\\")
    tex.append(f"{escape_char(what_sample)} read数 & {escape_char(mapped_line_num)} \\\\")
    tex.append("\\end{tabular}")
    return "%\n".join(tex)


def nanoplot_infos(nanoplot_qc_summary, nanoplot_raw_pic, nanoplot_clean_pic):
    tex = []
    tex.append("\\section{数据统计}")
    col_formats = (r">{\centering\arraybackslash}m{0.25\paperwidth}"
                   r">{\centering\arraybackslash}m{0.25\paperwidth}"
                   r">{\centering\arraybackslash}m{0.25\paperwidth}")
    tex.append(f"\\begin{{longtable}}{{{col_formats}}}")
    title = "\\hline 统计项 & 过滤前 & 过滤后 \\\\ \\hline"
    tex.append(title + "\\endfirsthead")
    tex.append(title + "\\endhead")
    tex.append("\\hline \\endfoot\n\\hline \\endlastfoot")
    df = pd.read_csv(nanoplot_qc_summary, header=None, sep="\t", skiprows=1)
    tex.append(df_to_table_tex(df))
    tex.append("\\end{longtable}\n")

    tex.append("\\begin{figure}[htbp]")
    tex.append("\\centering")
    tex.append("\\subfigure[过滤前]{\\includegraphics[width=0.5\\textwidth]{" + nanoplot_raw_pic + "}}")
    tex.append("\\subfigure[过滤后]{\\includegraphics[width=0.5\\textwidth]{" + nanoplot_clean_pic + "}}")
    tex.append("\\end{figure}")

    return "%\n".join(tex)


def coverage(cov_fig):
    tex = []
    tex.append("\\section{覆盖度}")
    tex.append("{\\raggedright \includegraphics[width=1.0\\textwidth]{" + cov_fig + "}}")
    return "%\n".join(tex)


def consensus(consensus_fa):
    tex = []
    tex.append("\\section{一致性序列}")
    tex.append(consensus_fa)
    return "%\n".join(tex)


def handle(unannotated_vcf):
    import vcf
    dir_name, file_name = path.split(unannotated_vcf)
    outfile = path.join(dir_name, f".unannotated.variants.info")
    with open(outfile, "w") as outf:
        variants = vcf.Reader(filename=unannotated_vcf)
        # the vcf may be annotated by medaka tools or called by longshot
        if "medaka_version" in dict(variants.metadata):
            for variant in variants:
                line = [str(variant.CHROM),
                        str(variant.POS),
                        f"{variant.REF}:{sum(variant.INFO['SR'][:2])}",
                        f"{variant.ALT[0]}:{sum(variant.INFO['SR'][-2:])}",
                        f"{variant.INFO['DP']}"]
                outf.write("\t".join(line) + "\n")
        else:
            for variant in variants:
                line = [str(variant.CHROM),
                        str(variant.POS),
                        f"{variant.REF}:{variant.INFO['AC'][0]}",
                        f"{variant.ALT[0]}:{variant.INFO['AC'][1]}",
                        f"{variant.INFO['DP']}"]
                outf.write("\t".join(line) + '\n')
    return outfile


def df_to_table_tex(df):
    lines = []
    for row in df.itertuples(index=False):
        new_row = map(lambda x: escape_char(x), row)
        lines.append(" & ".join(new_row) + "\\\\")
    return "%\n".join(lines)


def variants_list(annotated_var, unannotated_vcf):
    if not (path.exists(annotated_var) or path.exists(unannotated_vcf)):
        raise Exception()
    ann_col_formats = (r">{\centering\arraybackslash}m{0.08\paperwidth}"
                       r">{\centering\arraybackslash}m{0.06\paperwidth}"
                       r">{\centering\arraybackslash}m{0.06\paperwidth}"
                       r">{\centering\arraybackslash}m{0.06\paperwidth}"
                       r">{\centering\arraybackslash}m{0.06\paperwidth}"
                       r">{\centering\arraybackslash}m{0.1\paperwidth}"
                       r">{\centering\arraybackslash}m{0.08\paperwidth}"
                       r">{\centering\arraybackslash}m{0.08\paperwidth}"
                       r">{\centering\arraybackslash}m{0.15\paperwidth}"
                       r">{\centering\arraybackslash}m{0.15\paperwidth}")

    unann_col_formats = (r">{\centering\arraybackslash}m{0.15\paperwidth}"
                         r">{\centering\arraybackslash}m{0.15\paperwidth}"
                         r">{\centering\arraybackslash}m{0.15\paperwidth}"
                         r">{\centering\arraybackslash}m{0.15\paperwidth}"
                         r">{\centering\arraybackslash}m{0.15\paperwidth}")
    tex = []
    tex.append("\\section{变异列表}")
    tex.append("{\\noindent{\\scriptsize{")
    if path.exists(annotated_var):
        tex.append(f"\\begin{{longtable}}{{{ann_col_formats}}}")
        title = "\\hline 染色体 & 变异位置 & 参考碱基 & 突变碱基 & 位点深度 & 变异类型 & 变异影响 & 基因名 & 核苷酸变化 & 氨基酸变化 \\\\ \\hline"
        tex.append(title + "\\endfirsthead")
        tex.append(title + "\\endhead")
        tex.append("\\hline \\endfoot\n\\hline \\endlastfoot")
        df = pd.read_csv(annotated_var, sep="\t", header=None)
        tex.append(df_to_table_tex(df))
    else:
        tex.append(f"\\begin{{longtable}}{{{unann_col_formats}}}")
        title = "\\hline 染色体 & 变异位置 & 参考碱基 & 突变碱基 & 位点深度 \\\\ \\hline"
        tex.append(title + "\\endfirsthead")
        tex.append(title + "\\endhead")
        tex.append("\\hline \\endfoot\n\\hline \\endlastfoot")
        df = pd.read_csv(handle(unannotated_vcf), sep="\t", header=None)
        tex.append(df_to_table_tex(df))
    tex.append("\\end{longtable}\n}}}")
    return "%\n".join(tex)


def lineage(lineage_report):
    if not path.exists(lineage_report):
        return "%\n"
    from subprocess import run, PIPE

    tex = []
    result = run("pangolin --version", shell=True, stdout=PIPE)
    pangolin_version = result.stdout.decode("utf-8").strip()
    with open(lineage_report, 'r') as infile:
        lines = infile.readlines()
    contents = lines[1].split(",")
    pangolin_lineage, who = contents[1], contents[4]
    who = who if who else "\\_"
    tex.append("\\section{分型}")
    tex.append("\\begin{tabular}[htbp]{")
    tex.append(">{\\centering\\arraybackslash}m{0.2\\textwidth}" * 4 + "}\\hline")
    tex.append("样本 & Pangolin分型 & WHO命名 & Pangolin版本 \\\\ \\hline")
    tex.append(f"{sample_name} & {pangolin_lineage} & {who} & {pangolin_version} \\\\ \\hline")
    tex.append("\\end{tabular}")
    return "%\n".join(tex)


def tex_report(sample_name, out_tex):
    tex_lines = []
    files = get_files(sample_name)
    stat_fp = files['stat_fp']
    cov_fig_fp = files['cov_fig_fp']
    annotated_var_fp = files['annotated_var_list']
    unannotated_vcf_fp = files['unannotated_vcf_fp']
    consensus_fp = files['consensus_fp']
    lineage_report_fp = files['lineage_report_fp']
    nanoplot_qc_summary = files['nanoplot_qc_summary']
    nanoplot_raw_pic = files['nanoplot_raw_pic']
    nanoplot_clean_pic = files['nanoplot_clean_pic']

    tex_lines.append(preamble)
    tex_lines.append("\\begin{document}")
    tex_lines.append("\\begin{center}")
    tex_lines.append("\\LARGE{样本检测报告}")
    tex_lines.append("\\end{center}")
    tex_lines.append(infomations(stat_fp))
    tex_lines.append(nanoplot_infos(nanoplot_qc_summary, nanoplot_raw_pic, nanoplot_clean_pic))
    tex_lines.append(coverage(cov_fig_fp))
    tex_lines.append(consensus(consensus_fp))
    tex_lines.append(variants_list(annotated_var_fp, unannotated_vcf_fp))
    tex_lines.append(lineage(lineage_report_fp))
    tex_lines.append("\\end{document}")
    with open(out_tex, "w") as outfile:
        outfile.write("%\n".join(tex_lines))


if __name__ == '__main__':
    tex_report(sample_name, report_tex)
