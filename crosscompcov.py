#!/usr/bin/env python

"""
Get coverage upstream of a promoter.

Pass in a BED file with features to check, and a collection of read depth files
to check for read coverage statistics.

    python crosscompcov.py --bed promoters.bed organism/*.depths.txt > promoter.coverage.p
"""

import argparse
import pathlib
import re
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd


class Locations:
    def __init__(self, *, bed_file, bin_size, max_size, **kwargs):
        # read in the BED locations
        self.locations = pd.read_csv(bed_file, sep="\t", header=None, skiprows = 2, dtype={1:np.float64})
        # set column headers
        self.locations.columns = ["Chr", "Start", "Stop", "Gene"]
        self.bin_size = bin_size
        self.max_size = max_size

    # TODO: does this method name make sense?
    def binned_average(self, overlap, cutoff):
        return np.average(np.asarray(overlap[0:cutoff]).reshape(-1, self.bin_size), axis=1)

    def coverage_stats(self, feature_file):
        # read coverage file
        coverage = pd.read_csv(feature_file, sep="\t", header=None, dtype={1:np.float64})
        coverage.columns = ["Chr", "Start", "Cov"]
        # output values
        signal = []
        position = []
        average = []
        # loop the BED locations
        for _, row in self.locations.iterrows():
            # find overlapping coverage
            after_start = coverage["Start"] >= row["Start"]
            before_stop = coverage["Start"] <= row["Stop"]
            same_chr = coverage["Chr"] == row["Chr"]
            overlap = coverage[after_start & before_stop & same_chr]
            # pull out the columns we need
            sig = [r["Cov"] for _, r in overlap.iterrows()]
            pos = [r["Start"] for _, r in overlap.iterrows()]
            # calculate length to either the max size,
            # or largest that will fit evenly in bins
            cutoff = self.cutoff(len(overlap))
            binned_signal = self.binned_average(sig, cutoff)
            signal.append(list(binned_signal))
            binned_position = self.binned_average(pos, cutoff)
            position.append(list(binned_position))
            try:
                average.append(sum(pos) / len(pos))
            except ZeroDivisionError:
                average.append(0)
        return signal, position, average

    def cutoff(self, overlap_length):
        if overlap_length > self.max_size:
            return self.max_size
        else:
            return (overlap_length // self.bin_size) * self.bin_size


def setup_argument_parser():
    parser = argparse.ArgumentParser(
        description="Get coverage upstream of a promoter.",
        add_help=True,
    )
    # input files
    parser.add_argument("--bed-file", type=pathlib.Path)
    parser.add_argument("features", nargs="+", type=pathlib.Path)
    # output files
    parser.add_argument("--outfile", default=sys.stdout, type=pathlib.Path)
    # parameters
    parser.add_argument("--bin-size", default=100, type=int)
    parser.add_argument("--max-size", default=5000, type=int)
    return parser


if __name__ == "__main__":
    parser = setup_argument_parser()
    args = parser.parse_args()
    locs = Locations(**vars(args))
    # leave off process count so it can scale off the CPU of the host
    with Pool() as p:
        signal, position, average = zip(*p.map(locs.coverage_stats, args.features))
    # collect the stats by name
    columns = {}
    for f, s, p, a in zip(args.features, signal, position, average):
        # TODO support regex group to extract name from filename?
        # this is just removing all suffixes from the path
        name = f.name[:-len("".join(f.suffixes))]
        columns.update({
            f"{name}_signal": s,
            f"{name}_position": p,
            f"{name}_average": a,
        })
    # add stats to dataframe
    out_df = locs.locations.assign(**columns)
    # re-sort the columns for output
    _, i = np.unique(out_df.columns, return_index=True)
    pd.DataFrame(out_df.iloc[:, i]).to_pickle(args.outfile)
