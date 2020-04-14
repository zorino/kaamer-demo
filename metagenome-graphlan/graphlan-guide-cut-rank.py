#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-03-21
# version: 	0.01

import sys
import os

# Main #
if __name__ == "__main__":

    usage = """
    graphlan-guide-cut-rank.py <graphlan-guide.txt> <rank>


rank:
    phylum
    class
    order
    family
    genus

"""
    rank_symbols = ["nr__", "p__", "c__", "o__", "f__", "g__"]

    ranks = {
        "unknown": "nr__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__"
    }

    if len(sys.argv) < 3:
        print(usage)
        sys.exit(1)

    rank = sys.argv[2]
    rank_symbol = ranks[rank]
    saved_ranks = {}

    with open(sys.argv[1]) as f:
        for l in f:
            keep_rank = False
            current_rank = []
            for _r in l.rstrip().split("."):
                current_rank.append(_r)
                if rank_symbol in _r:
                    keep_rank = True
                    break

            if keep_rank:
                saved_ranks[".".join(current_rank)] = 1

    for _r, v in saved_ranks.items():
        print(_r)
