# This code covers Task 3 from Part 1

# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import glob
import subprocess

# configuration
BNFINDER_SCRIPT = os.path.join("tool", "bnf")
INPUT_DIR = "Project_Data_Final_3"
OUTPUT_DIR = "Inference_Results_3"
SCORING_CRITERIA = ["MDL", "BDE"]
# =================================================

def run_inference():
    """
    Function responsible for the reconstruction of boolean networks from the datasets generated in previous Task.
    For every data file, it applies BNFinder2 using both MDL and BDE scoring criteria.
    """
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    data_files = glob.glob(os.path.join(INPUT_DIR, "data_*.txt"))

    # loop iterating through all data files
    for data_file_path in data_files:
        filename = os.path.basename(data_file_path)
        file_root = os.path.splitext(filename)[0]

        print("\n--- Processing File: {} ---".format(filename))

        # using both 'MDL' and 'BDE' scoring criteria
        for score in SCORING_CRITERIA:
            output_filename = "result_{}_{}.sif".format(file_root, score)
            output_path = os.path.join(OUTPUT_DIR, output_filename)

            # constructing the command line arguments for BNFinder2
            # explanation:
            # -e: data file path
            # -s: scoring criterion
            # -n: output filename for the network
            # -v: prints progress
            # -l: max number of parents (in project the limit set to 3)
            cmd = [
                "python", BNFINDER_SCRIPT,
                "-e", data_file_path,
                "-s", score,
                "-n", output_path,
                "-v",
                "-l", "3"
            ]

            # executing the command
            subprocess.call(cmd)

if __name__ == "__main__":
    run_inference()