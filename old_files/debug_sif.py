# -*- coding: utf-8 -*-
import os
import glob

# Ustawienia
RESULTS_DIR = "Inference_Results_old"
TRUTH_DIR = "Project_Data_Final"


def check_files():
    print "--- DIAGNOSTYKA PLIKOW SIF ---"

    # 1. Sprawdzmy jeden plik .sif
    sif_files = glob.glob(os.path.join(RESULTS_DIR, "*.sif"))

    if not sif_files:
        print "[!] Nie znaleziono zadnych plikow .sif w folderze '{}'!".format(RESULTS_DIR)
        print "    Czy na pewno uruchomiles step3_inference.py z flaga -n?"
    else:
        # Bierzemy pierwszy z brzegu
        sif_file = sif_files[0]
        print "\n[1] Przykladowy plik SIF:", os.path.basename(sif_file)
        print "TRESC PLIKU (pierwsze 10 linii):"
        print "-" * 30
        try:
            with open(sif_file, 'r') as f:
                lines = f.readlines()
                if not lines:
                    print "[PUSTY PLIK]"
                for i, line in enumerate(lines[:10]):
                    print "LINE {}: {}".format(i, line.strip())
        except Exception as e:
            print "Blad odczytu: {}".format(e)

    # 2. Dla pewnosci sprawdzmy tez definicje
    print "\n[2] Przykladowy plik Ground Truth:"
    truth_files = glob.glob(os.path.join(TRUTH_DIR, "network_def_*.txt"))
    if truth_files:
        t_file = truth_files[0]
        print "Plik:", os.path.basename(t_file)
        with open(t_file, 'r') as f:
            for i, line in enumerate(f):
                if i < 3: print "LINE {}: {}".format(i, line.strip())


if __name__ == "__main__":
    check_files()