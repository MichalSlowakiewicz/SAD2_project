# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import csv
import random

# --- KONFIGURACJA ---
DATA_DIR = "Project_Data_Part2"
BASE_FILENAME = "data_neural_precursor_combined"
TRUTH_FILENAME = "true_network_structure.txt"
STEPS_list = [100, 200, 500, 1000]  # Dlugosci pojedynczej trajektorii
NUM_TRIALS = 10  # ILE ROZNYCH PROB DOKLEIC DO JEDNEGO PLIKU (np. 10 startow po 100 krokow)

if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)


def generate_clean_project():
    print("--- GENEROWANIE DANYCH (Wiele trajektorii w jednym pliku) ---")

    node_map = {
        1: "v_Neural_Precursor",
        2: "v_Progenitor",
        3: "v_elavl3_HuC",
        4: "v_her6",
        5: "v_miR_9",
        6: "v_zic5"
    }
    name_to_x = {v: "x{}".format(k) for k, v in node_map.items()}

    # Petla po dlugosci pojedynczej trajektorii
    for STEPS in STEPS_list:

        # Przygotowujemy slownik na dane z WIELU prob naraz
        # Na poczatku pusty
        combined_history = {name_to_x[name]: [] for name in node_map.values()}

        # Petla, ktora wykonuje X niezaleznych symulacji i skleja je
        for trial in range(NUM_TRIALS):

            # 1. Reset stanu (nowy losowy start)
            curr_state = {name: random.choice([0, 1]) for name in node_map.values()}

            # Dopisz stan poczatkowy tej proby do historii
            for name, val in curr_state.items():
                x_id = name_to_x[name]
                combined_history[x_id].append(val)

            # 2. Symulacja jednej sciezki
            for t in range(STEPS):
                s = curr_state
                n = {}

                # --- LOGIKA ---
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

                curr_state = n

                # Dopisz krok do historii
                for name, val in curr_state.items():
                    x_id = name_to_x[name]
                    combined_history[x_id].append(val)

        # --- ZAPIS (Jeden duzy plik po zakonczeniu wszystkich prob dla danej dlugosci) ---
        # Nazwa pliku np: data_combined_100.txt (co oznacza 10x100 krokow)
        DATA_FILENAME = "{}_{}.txt".format(BASE_FILENAME, STEPS)
        out_data = os.path.join(DATA_DIR, DATA_FILENAME)

        with open(out_data, 'wb') as f:
            writer = csv.writer(f, delimiter='\t')

            # Naglowek musi byc ciagly: 0, 1, 2 ... az do konca (NUM_TRIALS * (STEPS+1))
            total_len = len(combined_history["x1"])
            header = ["GeneName"] + [str(k) for k in range(total_len)]
            writer.writerow(header)

            for i in range(1, 7):
                x_key = "x{}".format(i)
                row = [x_key] + combined_history[x_key]
                writer.writerow(row)

        print("-> Wygenerowano zbiorczy plik: {} (Lacznie kolumn: {}, Sklada sie z {} trajektorii po {} krokow)".format(
            DATA_FILENAME, total_len, NUM_TRIALS, STEPS))

    # Zapis Prawdy
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


if __name__ == "__main__":
    generate_clean_project()