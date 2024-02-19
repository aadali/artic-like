from os import path
import math
from sys import argv

from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import pandas as pd


def xaxis_major_formater(x, pos=None):
    if x % 1000 == 0:
        return f"{int(x / 1000)}k"
    elif x % 100 == 0:
        return f"{x / 1000:.1f}k"


def plot_coverage(per_base_dp, figure_dir, sample_name):
    df = pd.read_csv(per_base_dp, sep="\t", header=None)
    df.columns = ['contig', 'position', 'depth']
    grouped = df.groupby("contig")
    for group in grouped.groups:
        sub_df = grouped.get_group(group)
        length = sub_df.shape[0]

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

        if length // 10000 > 10:  # For monkeypox whose length is 200k bp
            bin_width = 10000

        title = []
        for i in [10, 20, 50, 100, 500, 1000]:
            cov = sub_df.query("depth >= @i").shape[0] / length
            title.append(f">={i}X: {cov:.2%}")

        plot_title = f"{sample_name}-{group}\n{title[0]:^20}{title[1]:^20}{title[2]:^20}\n{title[3]:^20}{title[4]:^20}{title[5]:^20}"
        # plot_title = f"{sample_name}-{group}\n{title[0]:^20}{title[1]:^20}\n{title[2]:^20}{title[3]:^20}"
        fig, ax = plt.subplots(figsize=(12, 6))
        positions = sub_df['position']
        depth = sub_df['depth']

        depth2 = depth
        ax.bar(positions, depth2, width=1)
        ax.set_xlabel('Position')
        ax.set_ylabel("Depth")
        ax.set_title(plot_title)
        ax.xaxis.set_major_locator(MultipleLocator(bin_width))
        # ax.yaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_major_formatter(xaxis_major_formater)
        ax.set_xlim(-length * 0.01, length * 1.01)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(axis="x", labelsize=8)
        # plt.savefig(path.join(figure_dir, f"{sample_name}-{group}.raw.coverage.pdf"))
        plt.savefig(path.join(figure_dir, f"{sample_name}-{group}.raw.coverage.png"))

        fig, ax = plt.subplots(figsize=(12, 6))
        depth2 = [dp if dp <= 1050 else 1050 + (math.log(dp - 1050, 2) * 5) for dp in depth]
        ax.bar(positions, depth2, width=1)
        ax.set_xlabel('Position')
        ax.set_ylabel("Depth")
        ax.set_title(plot_title)
        ax.xaxis.set_major_locator(MultipleLocator(bin_width))
        # ax.yaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_major_formatter(xaxis_major_formater)
        ax.set_xlim(-length * 0.01, length * 1.01)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(axis="x", labelsize=8)
        # plt.savefig(path.join(figure_dir, f"{sample_name}-{group}.modified.coverage.pdf"))
        plt.savefig(path.join(figure_dir, f"{sample_name}-{group}.modified.coverage.png"))


if __name__ == '__main__':
    usage = f"usage: {path.basename(__file__)} <per_base_dp.bed> <fig_dir> <sample_name>"
    if len(argv) != 4:
        raise Exception(usage + "\n")
    per_base_dp = argv[1]
    figure_dir = argv[2]
    sample = argv[3]
    plot_coverage(per_base_dp, figure_dir, sample)
