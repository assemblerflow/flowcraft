#!/usr/bin/env python3

import sys
import json

core_file, f1, f2 = sys.argv[1:4]

try:
    sample_id = sys.argv[4]
except IndexError:
    sample_id = None


def get_core_genes(core_file):

    with open(core_file) as fh:
        core_genes = [x.strip() for x in fh.readlines()[1:]
                      if x.strip() != ""]

    return core_genes


def filter_core_genes(locus_array, info_array, core_genes):

    core_array = []

    for gene, info in zip(*[info_array, locus_array]):
        if gene in core_genes:
            core_array.append(info)

    return core_array


def assess_quality(core_array, core_genes):

    # Get the total number of missing loci. The sum/map approach aggretates
    # the sum of all possible missing loci symbols.
    missing_loci = ["LNF", "PLOT3", "PLOT5", "NIPH", "ALM", "ASM"]
    locus_not_found = sum(map(core_array.count, missing_loci))

    perc = float(locus_not_found) / float(len(core_genes))

    # Fail sample with higher than 2% missing loci
    with open(".status", "w") as fh:
        if perc > 0.02:
            status = "fail"
        elif perc > 0.003:
            status = "warning"
        else:
            status = "pass"

        fh.write(status)

    return status, perc


def get_table_data(data_obj, sample_id=None):

    header_map = dict((p, h) for p, h in enumerate(data_obj["header"]))
    table_data = []

    for sample, data in data_obj.items():

        if sample == "header":
            continue

        cur_data = []
        for pos, d in enumerate(data):
            cur_data.append({
                "header": header_map[pos],
                "value": d,
                "table": "chewbbaca"
            })

        table_data.append({
            "sample": sample_id if sample_id else sample,
            "data": cur_data
        })

    return table_data


def main():
    core_genes = get_core_genes(core_file)

    with open(f1) as f1h, open(f2) as f2h:

        j1 = json.load(f1h)
        j2 = json.load(f2h)

        sample_info = [(k, v) for k, v in j1.items() if "header" not in k]
        current_array = j1["header"]
        status_info = []
        for sample, info in sample_info:

            sample_name = sample_id if sample_id else sample

            core_results = filter_core_genes(info, current_array, core_genes)
            status, perc = assess_quality(core_results, core_genes)
            status_info.append({
                "sample": sample_name,
                "status": status,
                "lnfPercentage": perc
            })

        table_data = get_table_data(j2, sample_name)
        res = {"cagao": [j1, j2], "status": status_info,
               "tableRow": table_data}

        with open(".report.json", "w") as fh:
            fh.write(json.dumps(res, separators=(",", ":")))


main()
