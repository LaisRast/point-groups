#!/usr/bin/env python 
# coding: utf-8

## Generates the following two sequences to be used in OEIS:
## * number of 4-dimensional point groups of order n, and
## * number of 4-dimensional chiral point groups of order n.
## The precomputed catalogs are generated using the main 
## function in `counting_groups.sage` with `orderbound = 10000`.

from pathlib import Path
import json

## output directory
output_dir = "oeis/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

## generate data
for typ in ["all", "chiral"]:
    with open(f"precomputed/order_catalog_{typ}_orderbound_10000.json") as f:
        data = json.load(f)

    ## print DATA field
    print(f"{typ}: ")
    seq = f"{len(data['1'])}" if '1' in data else "0"
    for n in range(2, 66):
        new_seq = f"{seq}, {len(data[str(n)])}" if str(n) in data else f"{seq}, 0"
        # if len(new_seq) > 260: break
        seq = new_seq
    print(f"{typ.upper()} DATA (n_terms: {n}, length: {len(seq)})")
    print(seq)
    print()

    ## write b-file up to n=10000
    with open(f"{output_dir}/b-file_{typ}.txt", "w") as g:
        for n in range(1, 10001):
            _temp = g.write(f"{n} {len(data[str(n)]) if str(n) in data else 0}\n")
