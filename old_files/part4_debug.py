import os
import glob

# Ustawiamy sciezki (tak jak w glownym skrypcie)
TRUTH_DIR = "Project_Data_Final"
RESULTS_DIR = "Inference_Results_old"


def peek_files():
    print "--- DIAGNOSTYKA PLIKOW ---"

    # 1. Sprawdzamy jeden plik z definicja (Ground Truth)
    truth_files = glob.glob(os.path.join(TRUTH_DIR, "network_def_*.txt"))
    if truth_files:
        t_file = truth_files[0]
        print "\n[1] Plik Definicji (Ground Truth):", os.path.basename(t_file)
        print "Tresc (pierwsze 5 linii):"
        with open(t_file, 'r') as f:
            for i, line in enumerate(f):
                if i < 5: print "  LINE {}: {}".format(i, line.strip())
    else:
        print "\n[!] Nie znaleziono plikow network_def_*.txt!"

    # 2. Sprawdzamy jeden plik z wynikiem (BNFinder Result)
    result_files = glob.glob(os.path.join(RESULTS_DIR, "result_*.txt"))
    if result_files:
        r_file = result_files[0]
        print "\n[2] Plik Wynikowy (BNFinder):", os.path.basename(r_file)
        print "Tresc (pierwsze 5 linii):"
        with open(r_file, 'r') as f:
            for i, line in enumerate(f):
                if i < 5: print "  LINE {}: {}".format(i, line.strip())
    else:
        print "\n[!] Nie znaleziono plikow result_*.txt!"


if __name__ == "__main__":
    peek_files()