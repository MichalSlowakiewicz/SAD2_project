# -*- coding: utf-8 -*-
from __future__ import print_function
import networkx as nx
import boolean
import random
import itertools
import pandas as pd
import os
import glob
import subprocess
import re
import csv

# --- Configuration ---
DATA_DIR = "part1_generated_networks"
RESULTS_DIR = "part1_inference"
REPORT_FILE = "Final_Report_Part1.csv"
BNFINDER_SCRIPT = os.path.join("tool", "bnf")
SCORING_CRITERIA = ["MDL", "BDE"]

class BN(object):
    __bool_algebra = boolean.BooleanAlgebra()

    def __init__(self, list_of_nodes, list_of_functions):
        self.num_nodes = len(list_of_nodes)
        self.node_names = list_of_nodes
        self.list_of_nodes = []
        for node_name in list_of_nodes:
            node = self.__bool_algebra.Symbol(node_name)
            self.list_of_nodes.append(node)
        self.functions = []
        for fun in list_of_functions:
            self.functions.append(self.__bool_algebra.parse(fun).simplify())

    def __int_to_state(self, x):
        binary_str = format(x, '0' + str(self.num_nodes) + 'b')
        state = [int(char) for char in binary_str]
        return tuple(state)

    def evaluate_node(self, node_idx, state_tuple):
        substitution = {}
        for i in range(self.num_nodes):
            if state_tuple[i] == 1:
                substitution[self.list_of_nodes[i]] = self.__bool_algebra.TRUE
            else:
                substitution[self.list_of_nodes[i]] = self.__bool_algebra.FALSE
        result = self.functions[node_idx].subs(substitution).simplify()
        return 1 if result == self.__bool_algebra.TRUE else 0

    def get_neighbor_states(self, state):
        neighbors = set()
        for i in range(self.num_nodes):
            new_val = self.evaluate_node(i, state)
            next_state_list = list(state)
            if next_state_list[i] != new_val:
                next_state_list[i] = new_val
                neighbors.add(tuple(next_state_list))
            else:
                neighbors.add(state)
        return neighbors

    def generate_state_transition_system(self):
        G = nx.DiGraph()
        for i in range(2 ** self.num_nodes):
            current_state = self.__int_to_state(i)
            G.add_node(current_state)
            neighbors = self.get_neighbor_states(current_state)
            for neighbor in neighbors:
                G.add_edge(current_state, neighbor)
        return G

    def get_next_state_synchronous(self, state):
        next_state_list = []
        for i in range(self.num_nodes):
            next_state_list.append(self.evaluate_node(i, state))
        return tuple(next_state_list)

    def simulate_trajectory(self, start_state, steps, mode='sync', sampling_freq=1):
        trajectory = []
        current_state = start_state
        trajectory.append(current_state)
        for t in range(1, steps * sampling_freq + 1):
            if mode == 'sync':
                current_state = self.get_next_state_synchronous(current_state)
            else:
                node_to_update = random.randint(0, self.num_nodes - 1)
                new_val = self.evaluate_node(node_to_update, current_state)
                next_s = list(current_state)
                next_s[node_to_update] = new_val
                current_state = tuple(next_s)
            if t % sampling_freq == 0:
                trajectory.append(current_state)
        return trajectory

    def analyze_attractor_structure(self, start_state, max_steps=10000):
        history = []
        seen = {}
        curr = start_state
        for i in range(max_steps):
            if curr in seen:
                start_index = seen[curr]
                transient_part = history[:start_index]
                attractor_part = history[start_index:]
                return transient_part, attractor_part
            seen[curr] = i
            history.append(curr)
            curr = self.get_next_state_synchronous(curr)
        return history, []

    def simulate_trajectory_async(self, start_state, steps, sampling_freq=1):
        trajectory = [start_state]
        current_state = start_state
        for t in range(1, steps * sampling_freq + 1):
            node_to_update = random.randint(0, self.num_nodes - 1)
            new_val = self.evaluate_node(node_to_update, current_state)
            next_s = list(current_state)
            next_s[node_to_update] = new_val
            current_state = tuple(next_s)
            if t % sampling_freq == 0:
                trajectory.append(current_state)
        return trajectory

def generate_random_bn_strings(num_nodes, max_parents=3):
    node_names = ['x{}'.format(i + 1) for i in range(num_nodes)]
    functions = []
    for i in range(num_nodes):
        k = random.randint(1, min(max_parents, num_nodes))
        parents = random.sample(node_names, k)
        parent_combinations = list(itertools.product([0, 1], repeat=k))
        true_combinations = [c for c in parent_combinations if random.random() > 0.5]

        if not true_combinations:
            functions.append("0")
        elif len(true_combinations) == len(parent_combinations):
            functions.append("1")
        else:
            clauses = []
            for combo in true_combinations:
                literals = []
                for parent_idx, val in enumerate(combo):
                    p_name = parents[parent_idx]
                    literals.append(p_name if val == 1 else "~{}".format(p_name))
                clauses.append("(" + " & ".join(literals) + ")")
            full_func = " | ".join(clauses)
            functions.append(full_func)
    return node_names, functions

def save_to_file(data, colnames, filename):
    if not data: return
    df = pd.DataFrame(data, columns=colnames)
    df_transposed = df.T
    df_transposed.index.name = "GeneName"
    df_transposed.to_csv(filename, sep='\t', header=True, index=True)

def task_generate_data():
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)
    
    sizes_to_test = [5, 7, 10, 16]
    NUM_STARTS = 300

    print("--- Starting Data Generation ---")
    for N in sizes_to_test:
        print("Processing Network Size: {}".format(N))
        ASYNC_STEPS = N * 5
        names, funcs = generate_random_bn_strings(num_nodes=N, max_parents=3)
        bn = BN(names, funcs)

        with open("{}/network_def_{}.txt".format(DATA_DIR, N), "w") as f:
            for n_name, f_func in zip(names, funcs):
                f.write("{} = {}\n".format(n_name, f_func))

        # Single Trajectory
        start_node_single = tuple([random.randint(0, 1) for _ in range(N)])
        
        # Sync
        s_trans, s_attr = bn.analyze_attractor_structure(start_node_single)
        if s_attr:
            s_full = s_trans + (s_attr * 5)
            save_to_file(s_full, names, "{}/data_size{}_sync_SINGLE_full.txt".format(DATA_DIR, N))
            save_to_file(s_trans + s_attr[:3], names, "{}/data_size{}_sync_SINGLE_trans_heavy.txt".format(DATA_DIR, N))
            save_to_file(s_trans[-3:] + s_attr, names, "{}/data_size{}_sync_SINGLE_attr_heavy.txt".format(DATA_DIR, N))
            save_to_file(s_full[::3], names, "{}/data_size{}_sync_SINGLE_sampled_f3.txt".format(DATA_DIR, N))

        # Async
        s_raw_async = bn.simulate_trajectory_async(start_node_single, steps=ASYNC_STEPS)
        save_to_file(s_raw_async, names, "{}/data_size{}_async_SINGLE_full.txt".format(DATA_DIR, N))
        save_to_file(s_raw_async[:int(ASYNC_STEPS*0.3)], names, "{}/data_size{}_async_SINGLE_trans_heavy.txt".format(DATA_DIR, N))
        save_to_file(s_raw_async[-int(ASYNC_STEPS*0.3):], names, "{}/data_size{}_async_SINGLE_attr_heavy.txt".format(DATA_DIR, N))
        save_to_file(s_raw_async[::3], names, "{}/data_size{}_async_SINGLE_sampled_f3.txt".format(DATA_DIR, N))

        # Multiple Trajectories
        sync_full, sync_trans_heavy, sync_attr_heavy, sync_sampled_f3 = [], [], [], []
        current_start = tuple([random.randint(0, 1) for _ in range(N)])

        for _ in range(NUM_STARTS):
            trans, attr = bn.analyze_attractor_structure(current_start)
            if not attr:
                current_start = tuple([random.randint(0, 1) for _ in range(N)])
                continue
            full_traj = trans + (attr * 5)
            sync_full.extend(full_traj)
            sync_trans_heavy.extend(trans + attr[:3])
            sync_attr_heavy.extend(trans[-3:] + attr)
            sync_sampled_f3.extend(full_traj[::3])
            
            last_state = list(full_traj[-1])
            idx_to_flip = random.randint(0, N - 1)
            last_state[idx_to_flip] = 1 - last_state[idx_to_flip]
            current_start = tuple(last_state)

        save_to_file(sync_full, names, "{}/data_size{}_sync_full.txt".format(DATA_DIR, N))
        save_to_file(sync_trans_heavy, names, "{}/data_size{}_sync_transient_heavy.txt".format(DATA_DIR, N))
        save_to_file(sync_attr_heavy, names, "{}/data_size{}_sync_attractor_heavy.txt".format(DATA_DIR, N))
        save_to_file(sync_sampled_f3, names, "{}/data_size{}_sync_sampled_f3.txt".format(DATA_DIR, N))

        async_full, async_trans_heavy, async_attr_heavy, async_sampled_f3 = [], [], [], []
        
        for _ in range(NUM_STARTS):
            start_node = tuple([random.randint(0, 1) for _ in range(N)])
            raw_traj = bn.simulate_trajectory_async(start_node, steps=ASYNC_STEPS, sampling_freq=1)
            async_full.extend(raw_traj)
            async_sampled_f3.extend(raw_traj[::3])
            async_trans_heavy.extend(raw_traj[:int(ASYNC_STEPS*0.3)])
            async_attr_heavy.extend(raw_traj[-int(ASYNC_STEPS*0.3):])

        save_to_file(async_full, names, "{}/data_size{}_async_full.txt".format(DATA_DIR, N))
        save_to_file(async_trans_heavy, names, "{}/data_size{}_async_trans_heavy.txt".format(DATA_DIR, N))
        save_to_file(async_attr_heavy, names, "{}/data_size{}_async_attr_heavy.txt".format(DATA_DIR, N))
        save_to_file(async_sampled_f3, names, "{}/data_size{}_async_sampled_f3.txt".format(DATA_DIR, N))

def task_run_inference():
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    data_files = glob.glob(os.path.join(DATA_DIR, "data_*.txt"))
    print("\n--- Starting Inference ---")
    
    for data_file_path in data_files:
        filename = os.path.basename(data_file_path)
        file_root = os.path.splitext(filename)[0]
        print("Processing: {}".format(filename))

        for score in SCORING_CRITERIA:
            output_filename = "result_{}_{}.sif".format(file_root, score)
            output_path = os.path.join(RESULTS_DIR, output_filename)
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
    edges = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line: continue
            target, func = line.split("=", 1)
            target = target.strip()
            parents = set(re.findall(r'x\d+', func))
            if target in parents: parents.remove(target)
            edges[target] = parents
    return edges

def parse_sif_result(filepath):
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
    union_count = TP + FP + FN
    if union_count > 0:
        jaccard_index = float(TP) / union_count
        jaccard_distance = 1.0 - jaccard_index
    else:
        jaccard_distance = 0.0 if (len(truth_edges) == 0 and len(result_edges) == 0) else 1.0

    return TP, FP, FN, precision, recall, f1, shd, jaccard_distance

def task_run_evaluation():
    print("\n--- Starting Evaluation ---")
    with open(REPORT_FILE, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(
            ["Network Size", "Dataset Type", "Scoring", "TP", "FP", "FN", "Precision", "Recall", "F1_Score", "SHD", "Jaccard Distance"]
        )
        result_files = glob.glob(os.path.join(RESULTS_DIR, "result_*.sif"))
        result_files.sort()

        for r_file in result_files:
            filename = os.path.basename(r_file)
            core = filename.replace("result_data_", "").replace(".sif", "")
            parts = core.split("_")
            size_val = parts[0].replace("size", "")
            scoring = parts[-1]
            dataset_type = "_".join(parts[1:-1])

            truth_file = os.path.join(DATA_DIR, "network_def_{}.txt".format(size_val))
            truth_graph = parse_ground_truth(truth_file)
            result_graph = parse_sif_result(r_file)
            tp, fp, fn, prec, rec, f1, shd, jacc_dist = calculate_metrics(truth_graph, result_graph)

            writer.writerow([size_val, dataset_type, scoring, tp, fp, fn, "%.4f" % prec, "%.4f" % rec, "%.4f" % f1, int(shd), "%.4f" % jacc_dist])
            if f1 > 0:
                print(" -> [SUCCESS] {} | F1-Score: {:.2f}".format(filename, f1))
            else:
                print(" -> [NOTICE] {} | No edges correctly identified.".format(filename))

    print("\nEvaluation saved to: {}".format(REPORT_FILE))

if __name__ == "__main__":
    # Task 1 & 2: Generate Data
    task_generate_data()
    
    # Task 3: Inference
    task_run_inference()
    
    # Task 4: Evaluation
    task_run_evaluation()