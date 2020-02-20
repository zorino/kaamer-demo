#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-02-19
# version: 	0.01

import sys
import json


def load_taxa(eggNog_taxa_file):
    taxa = {}
    with open(eggNog_taxa_file) as f:
        for l in f:
            lA = l.split("\t")
            taxa[lA[0]] = lA[1]

    return taxa


def add_count(rank_counts, taxon, count):
    if taxon in rank_counts:
        rank_counts[taxon] += count
    else:
        rank_counts[taxon] = count


def process_taxa(res_file, taxa):

    rank_codes = ["sk", "p", "c", "o", "f", "g"]
    rank_names = [
        "superkingdom", "phylum", "class", "order", "family", "genus"
    ]
    rank_counts = {
        "sk": {
            "unknown": 0
        },
        "p": {
            "unknown": 0
        },
        "c": {
            "unknown": 0
        },
        "o": {
            "unknown": 0
        },
        "f": {
            "unknown": 0
        },
        "g": {
            "unknown": 0
        }
    }

    with open(res_file) as f:
        f.readline()  # skip header
        for l in f:
            lA = l.split("\t")
            taxon = taxa[lA[0]]
            count = int(lA[1])
            rank_seen = {}
            for t in taxon.split("."):
                __t = t.split("__")
                t_code = __t[0]
                t_name = "__".join(__t[1:])
                rank_seen[t_code] = 1
                if t_code in rank_codes:
                    add_count(rank_counts[t_code], t_name, count)
            for rk in rank_codes:
                if rk not in rank_seen:
                    rank_counts[rk]["unknown"] += count

    print(json.dumps(rank_counts, indent=2, sort_keys=False))

    res_file_prefix = res_file.replace(".tsv", "")
    for i, r in enumerate(rank_codes):
        output_file = res_file_prefix + "." + rank_names[i] + ".tsv"
        f = open(output_file, "w")
        total_count = sum(rank_counts[r].values())
        unknown = rank_counts[r]["unknown"]
        total_known = total_count - unknown

        sorted_counts = {
            k: v
            for k, v in sorted(
                rank_counts[r].items(), key=lambda item: item[1], reverse=True)
        }

        f.write("Feature\tCount\tRelAbun\tRelAbun_Known\n")
        f.write("%s\t%d\t%f\t%f\n" %
                ("unknown", unknown, float(unknown) / total_count, 0))
        for k, v in sorted_counts.items():
            if k != "unknown":
                f.write("%s\t%d\t%f\t%f\n" %
                        (k, v, float(v) / total_count, float(v) / total_known))
        f.close()


# Main #
if __name__ == "__main__":

    usage = """

metagenome-taxonomy.py <res.taxa.tsv> <uhgp-90_eggNOG.taxa.tsv>

# This will split the taxonomy results by taxonomic rank

# see https://gist.github.com/zorino/152a80cc343d15bda8810547a039e4b2 for uhgp-90_eggNOG.taxa.tsv

"""
    if len(sys.argv) < 3:
        print(usage)
        sys.exit(1)

    taxa = load_taxa(sys.argv[2])
    process_taxa(sys.argv[1], taxa)
