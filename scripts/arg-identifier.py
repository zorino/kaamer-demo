#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2020-01-22
# version: 	0.01

import sys


def parse_result(kaamer_res, identity_co, coverage_co):

    arg_hits = {}

    with open(kaamer_res) as f:

        header = f.readline().rstrip().split("\t")
        header_indices = {}

        for i, h in enumerate(header):
            header_indices[h] = i

        for l in f:
            fields = l.strip().split("\t")

            query_uid = fields[header_indices["QueryId"]]
            query_uid += "__"
            query_uid += fields[header_indices["QStart"]]
            query_uid += "__"
            query_uid += fields[header_indices["QEnd"]]

            if query_uid in arg_hits and arg_hits[query_uid][header_indices[
                    "Bitscore"]] > fields[header_indices["Bitscore"]]:
                continue

            if fields[header_indices["type"]] == "AMR" and fields[
                    header_indices["subtype"]] == "POINT":
                if float(fields[header_indices["%Identity"]]) >= 100:
                    arg_hits[query_uid] = fields

            if fields[header_indices["type"]] == "AMR" and fields[
                    header_indices["subtype"]] == "AMR":
                coverage = 100 * float(
                    fields[header_indices["AlnLength"]]) / float(
                        fields[header_indices["sequence_length"]])
                if float(fields[header_indices["%Identity"]]
                         ) >= identity_co and coverage >= coverage_co:
                    arg_hits[query_uid] = fields

    return header, arg_hits


# Main #
if __name__ == "__main__":

    usage = """
    python arg-identifier.py <kaamer ncbi_arg output> [%identity def=90] [%coverage def=90]

    POINT mutation require 100% identity

    Otherwise default % identity and % coverage are set at 90

    """

    if len(sys.argv) < 2:
        print(usage)
        exit(1)

    identity = 90
    coverage = 90
    if len(sys.argv) > 2:
        identity = float(sys.argv[2])
    if len(sys.argv) > 3:
        coverage = float(sys.argv[3])

    header, arg_hits = parse_result(sys.argv[1], identity, coverage)
    print("\t".join(header))
    for k in arg_hits:
        print("\t".join(arg_hits[k]))
