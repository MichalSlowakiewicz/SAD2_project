# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import csv
import random
import glob
import subprocess
import re

# --- Configuration ---
DATA_DIR = "part2_bio_network"
RESULTS_DIR = "part2_inference"
REPORT_FILE = "Final_Report_Part2.csv"
TRUTH_FILENAME = "true_network_structure.txt"
BNFINDER_SCRIPT = os.path.join("tool", "bnf")
SCORING_CRITERIA = ["MDL", "BDE"]

# Simulation Parameters
BASE_FILENAME = "data_neural_precursor_combined"
STEPS_LIST = [100, 200, 500, 1000]  # Length of single trajectory
NUM_TRIALS = 10  # Number of independent trials to combine per file

def task_generate_data():
    """
    Generates synthetic data for the Neural Precursor network using asynchronous updates.
    Combines multiple trials into single datasets for robustness.
    """
    print("--- Starting Data Generation (Part 2) ---")
    
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)

    # Node mapping for biological names to generic x-variables
    node_map = {
        1: "v_Neural_Precursor",
        2: "v_Progenitor",
        3: "v_elavl3_HuC",
        4: "v_her6",
        5: "v_miR_9",
        6: "v_zic5"
    }
    name_to_x = {v: "x{}".format(k) for k, v in node_map.items()}

    # Loop over different trajectory lengths
    for STEPS in STEPS_LIST:
        # Initialize dictionary to hold combined history for all trials
        combined_history = {name_to_x[name]: [] for name in node_map.values()}

        # Run multiple independent trials
        for trial in range(NUM_TRIALS):
            # 1. Random initialization
            curr_state = {name: random.choice([0, 1]) for name in node_map.values()}

            # Record initial state
            for name, val in curr_state.items():
                x_id = name_to_x[name]
                combined_history[x_id].append(val)

            # 2. Simulate trajectory (Asynchronous Update)
            for t in range(STEPS):
                s = curr_state
                n = {}

                # --- Biological Rules ---
                n['v_Neural_Precursor'] = s['v_elavl3_HuC']
                cond_prog = ((not s['v_her6'] and s['v_zic5']) or s['v_her6'])
                n['v_Progenitor'] = 1 if cond_prog else 0
                cond_huc = (not s['v_miR_9'] and not s['v_Progenitor'])
                n['v_elavl3_HuC'] = 1 if cond_huc else 0
                cond_her6 = (not s['v_miR_9'] and not s['v_Neural_Precursor'])
                n['v_her6'] = 1 if cond_her6 else 0
                cond_mir9 = (not s['v_her6'] and not s['v_Neural_Precursor'])
                n['v_miR_9'] = 1 if cond_mir9 else 0
                cond_zic5 = (not s['v_miR_9'] and not s['v_Neural_Precursor'])
                n['v_zic5'] = 1 if cond_zic5 else 0

                # Asynchronous update: pick one random node to update
                node_to_update = random.choice(list(node_map.values()))
                curr_state[node_to_update] = n[node_to_update]

                # Record state
                for name, val in curr_state.items():
                    x_id = name_to_x[name]
                    combined_history[x_id].append(val)

        # 3. Save combined dataset
        data_filename = "{}_{}.txt".format(BASE_FILENAME, STEPS)
        out_data = os.path.join(DATA_DIR, data_filename)

        with open(out_data, 'wb') as f:
            writer = csv.writer(f, delimiter='\t')

            # Create header: GeneName, 0, 1, 2, ...
            total_len = len(combined_history["x1"])
            header = ["GeneName"] + [str(k) for k in range(total_len)]
            writer.writerow(header)

            # Write rows for x1 to x6
            for i in range(1, 7):
                x_key = "x{}".format(i)
                row = [x_key] + combined_history[x_key]
                writer.writerow(row)

        print("-> Generated: {} (Total columns: {}, Trials: {}, Steps/Trial: {})".format(
            data_filename, total_len, NUM_TRIALS, STEPS))

    # Save Ground Truth Structure
    strict_rules = [
        "x1 = (x3)",
        "x2 = ((~x4 & x6) | x4)",
        "x3 = (~x5 & ~x2)",
        "x4 = (~x5 & ~x1)",
        "x5 = (~x4 & ~x1)",
        "x6 = (~x5 & ~x1)"
    ]
    out_truth = os.path.join(DATA_DIR, TRUTH_FILENAME)
    with open(out_truth, 'wb') as f:
        for rule in strict_rules:
            f.write(rule + "\n")
    print("-> Saved ground truth to: {}".format(TRUTH_FILENAME))


def task_run_inference():
    """
    Runs BNFinder on the generated datasets to reconstruct the network.
    """
    print("\n--- Starting Inference (Part 2) ---")
    
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    data_files = glob.glob(os.path.join(DATA_DIR, "data_*.txt"))
    
    if not data_files:
        print("Warning: No data files found in {}. Run generation first.".format(DATA_DIR))
        return

    for data_file_path in data_files:
        filename = os.path.basename(data_file_path)
        file_root = os.path.splitext(filename)[0]

        print("Processing: {}".format(filename))

        for score in SCORING_CRITERIA:
            output_filename = "result_{}_{}.sif".format(file_root, score)
            output_path = os.path.join(RESULTS_DIR, output_filename)

            # Command for BNFinder2
            # -l 3: limit parents to 3 (matching generation logic)
            cmd = [
                "python", BNFINDER_SCRIPT,
                "-e", data_file_path,
                "-s", score,
                "-n", output_path,
                "-v",
                "-l", "3"
            ]
            subprocess.call(cmd)


def parse_ground_truth(filepath):
    """
    Parses the ground truth file containing logical rules.
    Extracts parent nodes using regex.
    """
    edges = {}
    if not os.path.exists(filepath):
        print("Error: Truth file not found:", filepath)
        return edges

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line: continue

            target, func = line.split("=", 1)
            target = target.strip()

            # Extract variables like x1, x2, etc. ignoring operators
            parents = set(re.findall(r'x\d+', func))

            if target in parents: parents.remove(target)
            edges[target] = parents
    return edges


def parse_sif_result(filepath):
    """
    Parses BNFinder .sif output files.
    Format: Parent Interaction Child
    """
    edges = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2: continue

            parent = parts[0].strip()
            child = parts[-1].strip()

            if child not in edges: edges[child] = set()
            edges[child].add(parent)
    return edges


def calculate_metrics(truth_edges, result_edges):
    """
    Calculates TP, FP, FN, Precision, Recall, F1, SHD, and Jaccard Distance.
    """
    TP, FP, FN = 0, 0, 0
    all_nodes = set(truth_edges.keys()) | set(result_edges.keys())

    for node in all_nodes:
        true_parents = truth_edges.get(node, set())
        found_parents = result_edges.get(node, set())

        TP += len(true_parents.intersection(found_parents))
        FP += len(found_parents - true_parents)
        FN += len(true_parents - found_parents)

    precision = float(TP) / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = float(TP) / (TP + FN) if (TP + FN) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

    shd = FP + FN
    union = TP + FP + FN
    
    if union > 0:
        jaccard_index = float(TP) / union
        jaccard_dist = 1.0 - jaccard_index
    else:
        # If both are empty, distance is 0. If one is empty and other not, distance 1.
        jaccard_dist = 0.0 if (not truth_edges and not result_edges) else 1.0

    return TP, FP, FN, precision, recall, f1, shd, jaccard_dist


def task_run_evaluation():
    """
    Compares inference results against ground truth and saves a CSV report.
    """
    print("\n--- Starting Evaluation (Part 2) ---")

    truth_path = os.path.join(DATA_DIR, TRUTH_FILENAME)
    truth_graph = parse_ground_truth(truth_path)
    print("Loaded Ground Truth: {} nodes".format(len(truth_graph)))

    if not os.path.exists(RESULTS_DIR):
        print("Error: Results directory does not exist:", RESULTS_DIR)
        return

    result_files = glob.glob(os.path.join(RESULTS_DIR, "result_*.sif"))
    result_files.sort()

    if not result_files:
        print("No .sif result files found in:", RESULTS_DIR)
        return

    with open(REPORT_FILE, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(["Filename", "Traj_Length", "Scoring", "TP", "FP", "FN", "Precision", "Recall", "F1", "SHD", "Jaccard"])

        for r_file in result_files:
            filename = os.path.basename(r_file)

            # Determine scoring method from filename
            scoring = "MDL" if "MDL" in filename else "BDE"
            
            # Extract trajectory steps from filename (format: result_data_..._STEPS_SCORE.sif)
            # Example filename: result_data_neural_precursor_combined_100_MDL.sif
            # We split by '_' and grab the item before the score.
            parts = filename.replace(".sif", "").split("_")
            try:
                # Steps is the second to last element
                steps_val = int(parts[-2])
            except ValueError:
                steps_val = 0 # Fallback if naming convention differs

            result_graph = parse_sif_result(r_file)
            tp, fp, fn, prec, rec, f1, shd, jacc = calculate_metrics(truth_graph, result_graph)

            writer.writerow(
                [filename, steps_val, scoring, tp, fp, fn, "%.4f" % prec, "%.4f" % rec, "%.4f" % f1, shd, "%.4f" % jacc]
            )

            print(" -> {} (Steps: {}): F1 = {:.4f}".format(filename, steps_val, f1))

    print("\nReport saved to: {}".format(REPORT_FILE))


if __name__ == "__main__":
    # Task 1: Generate Data
    task_generate_data()
    
    # Task 2: Inference
    task_run_inference()
    
    # Task 3: Evaluation
    task_run_evaluation()