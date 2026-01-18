# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import glob
import re
import csv

# --- KONFIGURACJA DLA PART 2 ---
# Tutaj zmieniamy foldery na te, ktore stworzylismy w part2_clean_start.py
TRUTH_DIR = "Project_Data_Part2"
RESULTS_DIR = "Inference_Results_Part2"
REPORT_FILE = "Final_Report_Part2.csv"
TRUTH_FILENAME = "true_network_structure.txt"


def parse_ground_truth(filepath):
    """
    Ta funkcja jest SUPER - obsluguje format logiczny (x1 = ~x2 & x3)
    Dzieki re.findall wyciagnie same ID wezlow (x2, x3) ignorujac znaki specjalne.
    """
    edges = {}
    if not os.path.exists(filepath):
        print("BLAD: Nie znaleziono pliku prawdy:", filepath)
        return edges

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line: continue

            target, func = line.split("=", 1)
            target = target.strip()

            # To jest kluczowe: regex wyciaga 'x' z cyframi, ignorujac ~, &, |
            parents = set(re.findall(r'x\d+', func))

            if target in parents: parents.remove(target)
            edges[target] = parents
    return edges


def parse_sif_result(filepath):
    """
    Wczytuje wynik BNFindera.
    """
    edges = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2: continue

            # Format SIF: Rodzic +/- Dziecko
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
    union = TP + FP + FN
    jaccard = 1.0 - (float(TP) / union) if union > 0 else (0.0 if not truth_edges and not result_edges else 1.0)

    return TP, FP, FN, precision, recall, f1, shd, jaccard


def run_evaluation():
    print("--- EWALUACJA (PART 2) ---")

    # 1. Wczytujemy jeden, staly plik prawdy
    truth_path = os.path.join(TRUTH_DIR, TRUTH_FILENAME)
    truth_graph = parse_ground_truth(truth_path)
    print("Wczytano prawde: {} ({} wezlow)".format(truth_path, len(truth_graph)))

    # 2. Szukamy wynikow
    if not os.path.exists(RESULTS_DIR):
        print("BLAD: Folder wynikow nie istnieje:", RESULTS_DIR)
        return

    result_files = glob.glob(os.path.join(RESULTS_DIR, "result_*.sif"))
    result_files.sort()

    if not result_files:
        print("Nie znaleziono plikow .sif w:", RESULTS_DIR)
        return

    with open(REPORT_FILE, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(["Filename", "Scoring", "TP", "FP", "FN", "Precision", "Recall", "F1", "SHD", "Jaccard"])

        for r_file in result_files:
            filename = os.path.basename(r_file)

            # Proste wykrywanie metody (MDL/BDE) z nazwy pliku
            scoring = "MDL" if "MDL" in filename else "BDE"

            result_graph = parse_sif_result(r_file)
            tp, fp, fn, prec, rec, f1, shd, jacc = calculate_metrics(truth_graph, result_graph)

            writer.writerow(
                [filename, scoring, tp, fp, fn, "%.4f" % prec, "%.4f" % rec, "%.4f" % f1, shd, "%.4f" % jacc])

            print(" -> {} ({}): F1 = {:.4f}".format(filename, scoring, f1))

    print("\nRaport zapisany w: {}".format(REPORT_FILE))


if __name__ == "__main__":
    run_evaluation()