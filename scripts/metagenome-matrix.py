#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-02-17
# version: 	0.01

import sys
import os


def build_feature_matrix(ft, samples):

    ft_results = {}
    ft_results_count = {}

    for s in samples:
        with open("%s/res.%s.tsv" % (s, ft), "r") as f:
            f.readline()
            for l in f:
                lA = l.rstrip().split("\t")
                # print(lA)
                ft_key = lA[0].replace("unclassified ", "")

                if ft_key not in ft_results_count:
                    ft_results_count[ft_key] = 0

                if ft_key not in ft_results:
                    ft_results[ft_key] = {}

                ft_results[ft_key][s] = lA[3]
                ft_results_count[ft_key] += float(lA[3])

    sorted_fts = dict(
        sorted(ft_results_count.items(), key=lambda item: -item[1]))
    header = "Sample"
    for ft, v in sorted_fts.items():
        if ft == "unknown" or ft == "None":
            continue
        header += "\t%s" % ft

    print(header)

    for s in samples:
        output = os.path.basename(s)
        for ft, v in sorted_fts.items():
            if ft == "unknown" or ft == "None":
                continue
            if s in ft_results[ft]:
                output += "\t%s" % ft_results[ft][s]
            else:
                output += "\t0"
        print(output)


# Main #
if __name__ == "__main__":

    usage = """
python metagenome-matrix.py <ft> results-dir*

fts: go, ec, kegg_pathway, kegg_module, cog

fts (taxonomy):

  - taxa.superkingdom
  - taxa.phylum
  - taxa.class
  - taxa.order
  - taxa.family
  - taxa.genus

"""

if len(sys.argv) < 2:
    print(usage)
    sys.exit(1)

features = [
    "go",
    "ec",
    "kegg_pathway",
    "kegg_module",
    "taxa.superkingdom",
    "taxa.phylum",
    "taxa.class",
    "taxa.order",
    "taxa.family",
    "taxa.genus",
]

if sys.argv[1] not in features:
    print("\nUnrecognized feature!!\n")

build_feature_matrix(sys.argv[1], sys.argv[2:])
