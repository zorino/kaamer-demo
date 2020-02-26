#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-02-13
# version: 	0.01

import sys
import subprocess
import os
import gzip


def increment_counter(counter, keys):

    for k in keys:
        if k in counter:
            counter[k] += 1
        else:
            counter[k] = 1


def count_reads(fastq_file):

    nb_lines = -1
    if ".gz" in fastq_file:
        p1 = subprocess.Popen(["zcat", fastq_file], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["wc", "-l"],
                              stdin=p1.stdout,
                              stdout=subprocess.PIPE)
        out, err = p2.communicate()
        nb_lines = int(out)
    else:
        p1 = subprocess.Popen(["cat", fastq_file], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["wc", "-l"],
                              stdin=p1.stdout,
                              stdout=subprocess.PIPE)
        out, err = p2.communicate()
        nb_lines = int(out)

    return int(nb_lines) / 4


def process_result(fastq_file, res_file, output, min_identity, min_klength):

    total_reads = count_reads(fastq_file)
    total_reads_hit = 0
    go_counter = {}
    ec_counter = {}
    taxon_counter = {}
    kegg_pathway_counter = {}
    kegg_module_counter = {}
    cog_counter = {}

    if ".gz" in res_file:
        f = gzip.open(res_file, "rt")
    else:
        f = open(res_file, "r")

    f.readline()  # skip header

    last_id = ""
    for l in f:

        lA = l.rstrip('\n').split("\t")
        if last_id == lA[0]:
            continue

        last_id = lA[0]

        if float(lA[2]) < min_identity or float(lA[3]) < min_klength:
            continue

        total_reads_hit += 1
        taxon = lA[12]
        increment_counter(taxon_counter, [taxon])

        go = lA[13].split(",")
        increment_counter(go_counter, go)

        ec = lA[14].split(",")
        increment_counter(ec_counter, ec)

        kegg_pathway = lA[16].split(",")
        kegg_pathway = list(filter(lambda x: "map" not in x, kegg_pathway))
        increment_counter(kegg_pathway_counter, kegg_pathway)

        kegg_module = lA[17].split(",")
        increment_counter(kegg_module_counter, kegg_module)

        cog = list(lA[18])
        if len(cog) < 1:
            cog = [""]
        increment_counter(cog_counter, cog)

    print_results(
        {
            "total_reads": total_reads,
            "total_reads_hits": total_reads_hit,
            "go_counter": go_counter,
            "ec_counter": ec_counter,
            "taxon_counter": taxon_counter,
            "kegg_pathway_counter": kegg_pathway_counter,
            "kegg_module_counter": kegg_module_counter,
            "cog_counter": cog_counter
        }, output)


def print_results(results, output):

    print("Number of Reads: %d" % (results["total_reads"]))
    print("Number of Reads with Hits: %s (%f percent)" %
          (results["total_reads_hits"],
           (results["total_reads_hits"] * 100 / results["total_reads"])))

    if not os.path.exists(output):
        os.mkdir(output)

    output_summary = output + "/res.summary.tsv"

    # summary result
    with open(output_summary, "w") as f:
        f.write("Total Reads:\t%d\n" % (results["total_reads"]))
        f.write("Total Reads with Hits:\t%d\n" % (results["total_reads_hits"]))
        f.write("Percent Reads with Hits:\t%f\n" %
                (results["total_reads_hits"] * 100 / results["total_reads"]))
        f.write("GO unique hits:\t%d\n" % (len(results["go_counter"])))
        f.write("EC unique hits:\t%d\n" % (len(results["ec_counter"])))
        f.write("Taxon unique hits:\t%d\n" % (len(results["taxon_counter"])))
        f.write("Kegg_Pathway unique hits:\t%d\n" %
                (len(results["kegg_pathway_counter"])))
        f.write("Kegg_Module unique hits:\t%d\n" %
                (len(results["kegg_module_counter"])))
        f.write("Cog unique hits:\t%d\n" % (len(results["cog_counter"])))

    # go result
    output_go = output + "/res.go.tsv"
    sorted_go_counter = {
        k: v
        for k, v in sorted(results["go_counter"].items(),
                           key=lambda item: item[1],
                           reverse=True)
    }
    with open(output_go, "w") as f:
        f.write("GO Number\tReadsCount\tRelAbun\tRelAbun_Known\n")
        if "" in sorted_go_counter:
            unknown_count = sorted_go_counter[""]
        else:
            unknown_count = 0
        known_count = results["total_reads_hits"] - unknown_count
        for k, v in sorted_go_counter.items():
            relab_known = v / known_count
            if k == "":
                k = "None"
                relab_known = 0
            f.write("%s\t%d\t%f\t%f\n" %
                    (k, v, (v / results["total_reads_hits"]), relab_known))

    # ec result
    output_ec = output + "/res.ec.tsv"
    sorted_ec_counter = {
        k: v
        for k, v in sorted(results["ec_counter"].items(),
                           key=lambda item: item[1],
                           reverse=True)
    }
    with open(output_ec, "w") as f:
        f.write("EC Number\tReadsCount\tRelAbun\tRelAbun_Known\n")
        if "" in sorted_ec_counter:
            unknown_count = sorted_ec_counter[""]
        else:
            unknown_count = 0
        known_count = results["total_reads_hits"] - unknown_count
        for k, v in sorted_ec_counter.items():
            relab_known = v / known_count
            if k == "":
                k = "None"
                relab_known = 0
            f.write("%s\t%d\t%f\t%f\n" %
                    (k, v, (v / results["total_reads_hits"]), relab_known))

    # taxon result
    output_taxon = output + "/res.taxa.tsv"
    sorted_taxon_counter = {
        k: v
        for k, v in sorted(results["taxon_counter"].items(),
                           key=lambda item: item[1],
                           reverse=True)
    }
    with open(output_taxon, "w") as f:
        f.write("Taxon\tReadsCount\tRelAbun\tRelAbun_Known\n")
        if "" in sorted_taxon_counter:
            unknown_count = sorted_taxon_counter[""]
        else:
            unknown_count = 0
        known_count = results["total_reads_hits"] - unknown_count
        for k, v in sorted_taxon_counter.items():
            relab_known = v / known_count
            if k == "":
                k = "None"
                relab_known = 0
            f.write("%s\t%d\t%f\t%f\n" %
                    (k, v, (v / results["total_reads_hits"]), relab_known))

    # kegg_pathway result
    output_kegg_pathway = output + "/res.kegg_pathway.tsv"
    sorted_kegg_pathway_counter = {
        k: v
        for k, v in sorted(results["kegg_pathway_counter"].items(),
                           key=lambda item: item[1],
                           reverse=True)
    }
    with open(output_kegg_pathway, "w") as f:
        f.write("Kegg_Pathway\tReadsCount\tRelAbun\tRelAbun_Known\n")
        if "" in sorted_kegg_pathway_counter:
            unknown_count = sorted_kegg_pathway_counter[""]
        else:
            unknown_count = 0
        known_count = results["total_reads_hits"] - unknown_count
        for k, v in sorted_kegg_pathway_counter.items():
            relab_known = v / known_count
            if k == "":
                k = "None"
                relab_known = 0
            f.write("%s\t%d\t%f\t%f\n" %
                    (k, v, (v / results["total_reads_hits"]), relab_known))

    # kegg_module result
    output_kegg_module = output + "/res.kegg_module.tsv"
    sorted_kegg_module_counter = {
        k: v
        for k, v in sorted(results["kegg_module_counter"].items(),
                           key=lambda item: item[1],
                           reverse=True)
    }
    with open(output_kegg_module, "w") as f:
        f.write("Kegg_Module\tReadsCount\tRelAbun\tRelAbun_Known\n")
        if "" in sorted_kegg_module_counter:
            unknown_count = sorted_kegg_module_counter[""]
        else:
            unknown_count = 0
        known_count = results["total_reads_hits"] - unknown_count
        for k, v in sorted_kegg_module_counter.items():
            relab_known = v / known_count
            if k == "":
                k = "None"
                relab_known = 0
            f.write("%s\t%d\t%f\t%f\n" %
                    (k, v, (v / results["total_reads_hits"]), relab_known))

    # cog result
    output_cog = output + "/res.cog.tsv"
    sorted_cog_counter = {
        k: v
        for k, v in sorted(results["cog_counter"].items(),
                           key=lambda item: item[1],
                           reverse=True)
    }
    with open(output_cog, "w") as f:
        f.write("COG\tReadsCount\tRelAbun\tRelAbun_Known\n")
        if "" in sorted_cog_counter:
            unknown_count = sorted_cog_counter[""]
        else:
            unknown_count = 0
        known_count = results["total_reads_hits"] - unknown_count
        for k, v in sorted_cog_counter.items():
            relab_known = v / known_count
            if k == "":
                k = "None"
                relab_known = 0
            f.write("%s\t%d\t%f\t%f\n" %
                    (k, v, (v / results["total_reads_hits"]), relab_known))


# Main #
if __name__ == "__main__":

    usage = """
python metagenome-profiling.py <fastq> <kaamer-res> <output-dir> [%identity def=90] [min_query_klength def=10]
"""

    if len(sys.argv) < 4:
        print(usage)
        sys.exit(1)

    min_identity = 90
    min_klength = 10

    if len(sys.argv) > 4:
        min_identity = float(sys.argv[4])
    if len(sys.argv) > 5:
        min_klength = float(sys.argv[5])

    process_result(sys.argv[1], sys.argv[2], sys.argv[3], min_identity,
                   min_klength)
