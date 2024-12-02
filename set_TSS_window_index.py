#!/usr/bin/env python

"""
Process gene annotations to add buffer zones for promoters, on either side. If
input(s) are GFF3, convert to BED.

Using the gff2bed library, convert data from GFF3 to BED. List input file(s) as
arguments, and the outputs will be written to files with the BED extension:

    python set_TSS_window_index.py alpha.gff subdir/bravo.ext charlie.bed
    # will write: 
    #  alpha.bed
    #  alpha.promoter.bed
    #  subdir/bravo.promoter.bed 
    #  charlie.promoter.bed 
"""

import argparse
import csv
import itertools
import pathlib
import sys

import gff2bed


def convert(input_gff):
    return gff2bed.convert(gff2bed.parse(input_gff))


def extract_bed(path):
    if path.suffix in (".gff", ".gff3", ".gz"):
        return write_bed(path)
    return read_bed(path)


def filter_bed(bed, window):
    for row in bed:
        chrom, start, stop, gene, score, strand = row
        # choose window center based on whether positive or negative strand
        position = start if strand == "+" else stop
        # filter out any windows that go beyond start of strand
        if position - window > 0:
            yield chrom, position - window, position + window, gene


def process(filename, window):
    half = window // 2
    for path in filename:
        bed = extract_bed(path)
        write_promoter(path, filter_bed(bed, half))


def read_bed(path):
    with open(path, "rt") as f:
        reader = csv.reader(f, dialect="excel-tab")
        yield from reader


def setup_argument_parser():
    parser = argparse.ArgumentParser(
        description="Prepare annotations for promoter matching.",
        add_help=True,
    )
    parser.add_argument("filename", nargs="+", type=pathlib.Path)
    parser.add_argument("--window", default=5000, nargs=1, type=int)
    return parser


def write_bed(gff_path):
    with open(gff_path.parent / f"{gff_path.stem}.bed", "wt") as f:
        writer = csv.writer(f, dialect="excel-tab")
        for row in convert(gff_path):
            writer.writerow(row)
            yield row


def write_promoter(path, regions):
    with open(path.parent / f"{path.stem}.promoter.bed", "wt") as f:
        writer = csv.writer(f, dialect="excel-tab")
        writer.writerows(regions)


def main():
    parser = setup_argument_parser()
    args = parser.parse_args()
    process(**vars(args))


if __name__ == "__main__":
    main()
