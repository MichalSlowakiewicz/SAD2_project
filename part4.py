# This code covers Task 4 from Part 1

# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import glob
import re
import csv

# configuration
TRUTH_DIR = "Project_Data_Final"
RESULTS_DIR = "Inference_Results_old"
REPORT_FILE = "Final_Report_Part1.csv"

def parse_ground_truth(filepath):
    """
    Parses the original network definition file to extract true edges.
    Returns a dictionary: {child_node: set(parent_nodes)}.
    """
    edges = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line: continue

            target, func = line.split("=", 1)
            target = target.strip()

            # obtaining node names from the logical expression
            parents = set(re.findall(r'x\d+', func))

            # removing potential self-loops
            if target in parents: parents.remove(target)
            edges[target] = parents
    return edges


def parse_sif_result(filepath):
    """
    Parses BNFinder '.sif' output files.
    """
    edges = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2: continue

            # BNFinder2 has format: parent, interaction_type, child
            parent = parts[0].strip()
            child = parts[-1].strip()

            if child not in edges: edges[child] = set()
            edges[child].add(parent)
    return edges


def calculate_metrics(truth_edges, result_edges):
    """
    Calculates standard metrics: TP, FP, FN, Precision, Recall, F1,
    by comparing the results from Task 3 against the ground truth.
    """
    TP, FP, FN = 0, 0, 0
    all_nodes = set(truth_edges.keys()) | set(result_edges.keys())

    for node in all_nodes:
        true_parents = truth_edges.get(node, set())
        found_parents = result_edges.get(node, set())

        # true positives - correctly identified parents
        TP += len(true_parents.intersection(found_parents))
        # false positives - predicted parents that don't exist in truth
        FP += len(found_parents - true_parents)
        # false negatives - true parents that were missed
        FN += len(true_parents - found_parents)

    precision = float(TP) / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = float(TP) / (TP + FN) if (TP + FN) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

    return TP, FP, FN, precision, recall, f1


def run_evaluation():
    """
    Main evaluation pipeline - it compares the results from Task 3 with
    corresponding ground truths and saves final metrics to a '.csv' file.
    """
    with open(REPORT_FILE, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(
            ["Network Size", "Dataset Type", "Scoring", "TP", "FP", "FN", "Precision", "Recall", "F1_Score"]
        )

        result_files = glob.glob(os.path.join(RESULTS_DIR, "result_*.sif"))
        result_files.sort()

        for r_file in result_files:
            filename = os.path.basename(r_file)
            # parsing information from filename
            core = filename.replace("result_data_", "").replace(".sif", "")
            parts = core.split("_")

            size_val = parts[0].replace("size", "")
            scoring = parts[-1]
            dataset_type = "_".join(parts[1:-1])

            # locating corresponding ground truth file
            truth_file = os.path.join(TRUTH_DIR, "network_def_{}.txt".format(size_val))

            truth_graph = parse_ground_truth(truth_file)
            result_graph = parse_sif_result(r_file)

            # calculating metrics
            tp, fp, fn, prec, rec, f1 = calculate_metrics(truth_graph, result_graph)

            # saving results to '.csv'
            writer.writerow([size_val, dataset_type, scoring, tp, fp, fn, "%.4f" % prec, "%.4f" % rec, "%.4f" % f1])

            # displaying some information about results
            if f1 > 0:
                print("  -> [SUCCESS] {} | F1-Score: {:.2f}".format(filename, f1))
            else:
                print("  -> [NOTICE] {} | No edges correctly identified.".format(filename))

    # finish message
    print("\nEvaluation saved to: {}".format(REPORT_FILE))


if __name__ == "__main__":
    run_evaluation()