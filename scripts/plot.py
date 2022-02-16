from os import path
from sys import argv

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd


# per_base_dp = "/home/a/big/ycq/projects/artic-like/test001/mosdepth/test001.per-base.bed.gz"
# fig_path = "test001.coverage.jpg"


def coverage_ratio(df, depth):
    """
    calculate coverage percentage depend on the depth
    :param df:
    :param depth:
    :return:
    """
    sub_df = df.query("depth >= @depth")
    sub_df['count'] = df['end'] - df['start']
    length = df['end'].tolist()[-1]
    return sum(sub_df['count']) / length


def xaxis_major_formater(x, pos=None):
    if x % 1000 == 0:
        return f"{int(x / 1000)}k"
    elif x % 100 == 0:
        return f"{x / 1000:.1f}k"


def plot_coverage(per_base_dp, out_dir):
    df = pd.read_csv(per_base_dp, compression="gzip", sep="\t", header=None)
    df.columns = ['contig', 'start', 'end', 'depth']
    grouped = df.groupby("contig")
    for group in grouped.groups:
        sub_df = grouped.get_group(group)
        length = df['end'].tolist()[-1]

        # the different length contig has different bin_width to set ticker
        bins = length // 1000
        if bins > 10:
            bin_width = 1000
        elif 5 < bins <= 10:
            bin_width = 500
        elif 1 < bins <= 5:
            bin_width = 200
        else:
            bin_width = 50

        title = []
        for i in [1, 10, 100, 1000]:
            cov = coverage_ratio(sub_df, i)
            title.append(f">={i}X: {cov:.2%}")

        plot_title = f"{group}\n{title[0]:^20}{title[1]:^20}\n{title[2]:^20}{title[3]:^20}"
        fig, ax = plt.subplots(figsize=(12, 6))
        positions = []
        depth = []
        for row in sub_df.itertuples(index=False):
            positions.extend(list(range(row.start, row.end)))
            depth.extend((row.end - row.start) * [row.depth])

        depth2 = [dp if dp < 1200 else 1200 for dp in depth]
        ax.bar(positions, depth2, width=1)
        ax.set_xlabel('Position')
        ax.set_ylabel("Depth")
        ax.set_title(plot_title)
        ax.xaxis.set_major_locator(MultipleLocator(bin_width))
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.xaxis.set_major_formatter(xaxis_major_formater)
        ax.set_xlim(-length * 0.01, length * 1.01)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(axis="x", labelsize=8)
        plt.savefig(path.join(out_dir, f"{group}.coverage.pdf"))


if __name__ == '__main__':
    usage = f"usage: {path.basename(__file__)} <per_base_dp.bed> <out_fig_path>"
    if len(argv) != 3:
        raise Exception(usage + "\n")
    per_base_dp = argv[1]
    out_dir = argv[2]
    plot_coverage(per_base_dp, out_dir)
