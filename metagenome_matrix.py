#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-02-17
# version: 	0.01

import sys
import os


def build_feature_matrix(ft, output_dir, samples):

    ft_results = {}
    ft_results_count = {}
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

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

                if s not in ft_results[ft_key]:
                    ft_results[ft_key][s] = {}

                ft_results[ft_key][s]["count"] = lA[1]
                ft_results[ft_key][s]["relab"] = lA[3]
                ft_results_count[ft_key] += float(lA[3])

    sorted_fts = dict(
        sorted(ft_results_count.items(), key=lambda item: -item[1]))
    header = "Sample"
    for _ft, v in sorted_fts.items():
        if _ft == "unknown" or _ft == "None":
            continue
        header += "\t%s" % _ft

    count_file = open("%s/res.%s.count.tsv" % (output_dir, ft), "w")
    relab_file = open("%s/res.%s.relab.tsv" % (output_dir, ft), "w")

    count_file.write(header + "\n")
    relab_file.write(header + "\n")

    for s in samples:
        output_count = os.path.basename(s)
        output_relab = os.path.basename(s)
        for _ft, v in sorted_fts.items():
            if _ft == "unknown" or _ft == "None":
                continue
            if s in ft_results[_ft]:
                output_relab += "\t%s" % ft_results[_ft][s]["relab"]
                output_count += "\t%s" % ft_results[_ft][s]["count"]
            else:
                output_relab += "\t0"
                output_count += "\t0"

        output_relab += "\n"
        output_count += "\n"
        relab_file.write(output_relab)
        count_file.write(output_count)

    count_file.close()
    relab_file.close()


# Main #
if __name__ == "__main__":

    usage = """

python metagenome-matrix.py <ft> <output-dir> results-dir*

fts: go, ec, kegg_pathway, kegg_module, cog

fts (taxonomy):

  - taxa.superkingdom
  - taxa.phylum
  - taxa.class
  - taxa.order
  - taxa.family
  - taxa.genus


# Will produce three matrices :
 -- raw count matrix (without unknown)
 -- abundance matrix (without unknown)
 -- clr log-ratio transform matrix (Aitchison)

"""

if len(sys.argv) < 3:
    print(usage)
    sys.exit(1)

features = [
    "go", "ec", "kegg_pathway", "kegg_module", "taxa.superkingdom",
    "taxa.phylum", "taxa.class", "taxa.order", "taxa.family", "taxa.genus",
    "cog"
]

if sys.argv[1] not in features:
    print("\nUnrecognized feature!!\n")

build_feature_matrix(sys.argv[1], sys.argv[2], sys.argv[3:])
