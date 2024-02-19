import sys
import gzip
from os import path
from math import log

import numpy as np
import pandas as pd
from Bio import SeqIO
from nanomath import get_N50
import plotly.express as px

useage = f"python3 {path.basename(__file__)}  <input_fastq>  <out_fig_path.png> <out_stat_path.txt>"
color = "#4CB391"
colormap = "Greens"


def errs_tab(n):
    """Generate list of error rates for qualities less than equal than n."""
    return [10 ** (q / -10) for q in range(n + 1)]


def ave_qual(quals, qround=False, tab=errs_tab(128)):
    """Calculate average basecall quality of a read.
    Receive the integer quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    """
    if quals:
        mq = -10 * log(sum([tab[q] for q in quals]) / len(quals), 10)
        if qround:
            return round(mq)
        else:
            return mq
    # else:
    #     return None


def extract_from_fastq(fq):
    for record in SeqIO.parse(fq, "fastq"):
        yield ave_qual(record.letter_annotations["phred_quality"]), len(record)


def get_fastq_input(fastq):
    assert isinstance(fastq, str)
    if fastq.endswith(".gz"):
        input_fastq = gzip.open(fastq, "rt")
    else:
        input_fastq = open(fastq, 'r')

    return pd.DataFrame(
        data=[record for record in extract_from_fastq(input_fastq) if record],
        columns=["quals", "lengths"]
    ).dropna()


def scatter(x, y, fig_path, names=["Read length", "Average read quality"],
            title="Read lengths vs Average read quality plot using dots"):
    fig = px.scatter(x=x,
                     y=y,
                     marginal_x="histogram",
                     marginal_y="histogram",
                     range_x=[0, np.amax(x)],
                     range_y=[0, np.amax(y)])

    fig.update_traces(marker=dict(color=color))
    fig.update_yaxes(rangemode="tozero")
    fig.update_xaxes(rangemode="tozero")
    fig.update_layout(xaxis_title=names[0],
                      yaxis_title=names[1],
                      title=title,
                      title_x=0.5
                      )
    fig.write_image(fig_path)


def get_stat(df, outfile):
    assert isinstance(df, pd.DataFrame)
    content = []
    content.append(f"Metrics\tdataset")
    content.append(f"number_of_reads\t{df.shape[0]}")
    content.append(f"number_of_bases\t{sum(df['lengths'])}")
    content.append(f"median_read_length\t{np.median(df['lengths'])}")
    content.append(f"mean_read_length\t{np.mean(df['lengths']):.2f}")
    content.append(f"n50\t{get_N50(df['lengths'])}")
    content.append(f"mean_qual\t{np.mean(df['quals']):.2f}")
    content.append(f"median_qual\t{np.median(df['quals']):.2f}")
    longest = df.sort_values(by="lengths", ascending=False).head(5)
    for idx, row in enumerate(longest.itertuples(index=False)):
        content.append(f"longest_read_(with_Q):{idx + 1}\t{row.lengths} ({row.quals:.2f})")
    highest = df.sort_values(by="quals", ascending=True).head(5)
    for idx, row in enumerate(highest.itertuples(index=False)):
        content.append(f"higest_Q_read_(with_length):{idx + 1}\t{row.quals:.2f} ({row.lengths})")

    for qual in [5, 7, 10, 12, 15, 20]:
        sub_df = df.query("quals >= @qual").reset_index()
        sub_reads_num = sub_df.shape[0]
        sub_lengths = sum(sub_df['lengths'])
        proportion = f"{sub_lengths / sum(df['lengths']):.2%}"
        sub_lengths_mb = f"{sub_lengths / 1000000:.2f}Mb"
        content.append(f"Reads >Q{qual}\t{sub_reads_num} ({proportion}) {sub_lengths_mb}")

    with open(outfile, 'w') as outf:
        outf.write("\n".join(content) + "\n")


if __name__ == '__main__':
    if len(sys.argv) != 4:
        raise Exception(useage)
    input_fastq, fig_path, stat_path = sys.argv[1:]
    df = get_fastq_input(input_fastq)
    scatter(df['lengths'], df['quals'], fig_path)
    get_stat(df, stat_path)
