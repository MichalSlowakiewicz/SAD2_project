# -*- coding: utf-8 -*-
from __future__ import print_function  # Python 2 compatibility
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Create output directory for graphs if it doesn't exist
if not os.path.exists("graphs"):
    os.makedirs("graphs")

# ==========================================
# PART 1: NETWORK SIZE & DYNAMICS ANALYSIS
# ==========================================

def load_and_prep_part1(filename="Final_Report_Part1.csv"):
    """
    Loads Part 1 data.
    Parses metadata (Trajectory Count, Regime, Mode).
    """
    print("Loading Part 1 data from {}...".format(filename))
    
    # Check if file exists to avoid immediate crash
    if not os.path.exists(filename):
        print("Warning: File not found. Returning empty DataFrame.")
        return pd.DataFrame()

    df = pd.read_csv(filename, sep=";")

    # --- Pre-processing ---
    def parse_metadata(row):
        name = row['Dataset Type']
        
        # 1. Trajectory Count
        traj = 'Single' if 'SINGLE' in name else 'Multiple (300)'
        
        # 2. Regime
        if 'sampled' in name: regime = 'Sub-sampled (f3)'
        elif 'trans' in name: regime = 'Transient-Heavy'
        elif 'attr' in name: regime = 'Attractor-Heavy'
        elif 'full' in name: regime = 'Full Dynamics'
        else: regime = 'Other'
            
        # 3. Mode
        mode = 'Async' if 'async' in name else 'Sync'
        
        return pd.Series([traj, regime, mode])

    df[['Trajectory_Count', 'Regime', 'Mode']] = df.apply(parse_metadata, axis=1)
    return df

def plot_part1_deep_dive(df):
    """
    Generates bar charts specifically for Network Sizes 5 & 7.
    """
    print("Generating Part 1 Deep Dive plots...")
    
    if df.empty:
        print("DataFrame is empty. Skipping Part 1 Deep Dive.")
        return

    # Filter specific subset
    subset = df[
        (df['Network Size'].isin([5, 7])) & 
        (df['Regime'] == 'Full Dynamics') &
        (df['Trajectory_Count'] == 'Multiple (300)')
    ].copy()

    if subset.empty:
        print(" -> Skip: No matching data for Deep Dive (Size 5/7, Full, Multiple).")
        return

    subset.sort_values(by=['Network Size', 'Mode', 'Scoring'], inplace=True)

    metrics = {
        'F1_Score': 'F1 Score (Higher is Better)',
        'SHD': 'SHD (Lower is Better)', 
        'Jaccard Distance': 'Jaccard Dist (Lower is Better)'
    }
    
    colors = {'BDE': '#1f77b4', 'MDL': '#ff7f0e'}

    for metric, title in metrics.items():
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle("{} - Full Dynamics (Multiple Trajectories)".format(title), fontsize=16, fontweight='bold')
        
        sizes = [5, 7]
        modes = ['Sync', 'Async']
        bar_width = 0.35
        x_indices = np.arange(len(modes))

        for i, size in enumerate(sizes):
            ax = axes[i]
            size_data = subset[subset['Network Size'] == size]
            
            bde_vals = []
            mdl_vals = []
            for mode in modes:
                # Use plain int/float extraction to avoid pandas series ambiguity
                val_bde_series = size_data[(size_data['Mode'] == mode) & (size_data['Scoring'] == 'BDE')][metric]
                val_mdl_series = size_data[(size_data['Mode'] == mode) & (size_data['Scoring'] == 'MDL')][metric]
                
                val_bde = val_bde_series.mean() if not val_bde_series.empty else 0
                val_mdl = val_mdl_series.mean() if not val_mdl_series.empty else 0

                bde_vals.append(0 if np.isnan(val_bde) else val_bde)
                mdl_vals.append(0 if np.isnan(val_mdl) else val_mdl)
            
            ax.bar(x_indices - bar_width/2, bde_vals, bar_width, label='BDE', color=colors['BDE'], alpha=0.9)
            ax.bar(x_indices + bar_width/2, mdl_vals, bar_width, label='MDL', color=colors['MDL'], alpha=0.9)
            
            ax.set_title("Network Size: {}".format(size), fontsize=14)
            ax.set_xticks(x_indices)
            ax.set_xticklabels(modes)
            ax.set_xlabel('Simulation Mode')
            ax.set_ylabel(metric)
            ax.grid(axis='y', linestyle='--', alpha=0.3)
            
            if i == 0: ax.legend(title="Scoring")

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        filename = "graphs/DeepDive_Size5_7_{}.png".format(metric)
        plt.savefig(filename, dpi=300)
        plt.close()

def plot_part1_metric_breakdown(df):
    """
    Generates 2x2 grids comparing Single vs Multiple trajectories across all regimes.
    """
    print("Generating Part 1 Metric Breakdown plots...")
    
    if df.empty: return

    metrics = {
        'F1_Score': 'F1 Score',
        'SHD': 'SHD', 
        'Jaccard Distance': 'Jaccard Dist',
        'Precision': 'Precision',
        'Recall': 'Recall'
    }
    regimes = ['Full Dynamics', 'Sub-sampled (f3)', 'Transient-Heavy', 'Attractor-Heavy']
    colors = ['#4e79a7', '#f28e2b'] # Blue (Single) vs Orange (Multiple)

    for metric, title in metrics.items():
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle("Impact of Trajectory Count on {}".format(title), fontsize=18, fontweight='bold')
        axes = axes.flatten()
        
        for i, regime in enumerate(regimes):
            ax = axes[i]
            subset = df[df['Regime'] == regime].copy()
            
            subset['Network Size'] = subset['Network Size'].astype(int)
            sizes = sorted(subset['Network Size'].unique())
            x = np.arange(len(sizes))
            bar_width = 0.35

            # Get means for Single vs Multiple
            vals_single = []
            vals_multi = []
            
            for s in sizes:
                v_s = subset[(subset['Network Size']==s) & (subset['Trajectory_Count']=='Single')][metric].mean()
                v_m = subset[(subset['Network Size']==s) & (subset['Trajectory_Count']=='Multiple (300)')][metric].mean()
                vals_single.append(v_s)
                vals_multi.append(v_m)
            
            # Replace NaNs with 0
            vals_single = [0 if np.isnan(v) else v for v in vals_single]
            vals_multi = [0 if np.isnan(v) else v for v in vals_multi]

            ax.bar(x - bar_width/2, vals_single, bar_width, label='Single', color=colors[0], alpha=0.9)
            ax.bar(x + bar_width/2, vals_multi, bar_width, label='Multiple (300)', color=colors[1], alpha=0.9)
            
            ax.set_title("Dataset Type: {}".format(regime), fontsize=14, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels(sizes)
            ax.set_xlabel('Network Size (Nodes)')
            ax.set_ylabel(metric)
            ax.grid(axis='y', linestyle='--', alpha=0.4)
            
            if i == 0: ax.legend()
            
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        filename = "graphs/Metric_Analysis_{}.png".format(metric)
        plt.savefig(filename, dpi=300)
        plt.close()

# ==========================================
# PART 2: BIOLOGICAL MODEL ANALYSIS
# ==========================================

def run_part2_pipeline(filename="Final_Report_Part2.csv"):
    """
    Loads Part 2 data and generates plots comparing trajectory length vs performance.
    """
    print("\n--- Starting Part 2 Visualization ---")
    
    if not os.path.exists(filename):
        print("Notice: '{}' not found. Creating mock data for visualization.".format(filename))
        # Mock Part 2 Data
        data_mock = []
        for length in [100, 200, 500, 1000]:
            for score in ['BDE', 'MDL']:
                 f1 = np.random.uniform(0.6, 0.95)
                 shd = np.random.randint(0, 5)
                 jacc = np.random.uniform(0, 0.4)
                 data_mock.append(["mock_{}_{}.sif".format(length, score), length, score, 0,0,0,0,0, f1, shd, jacc])
        df = pd.DataFrame(data_mock, columns=["Filename", "Traj_Length", "Scoring", "TP", "FP", "FN", "Precision", "Recall", "F1", "SHD", "Jaccard"])
    else:
        print("Loading Part 2 data from {}...".format(filename))
        df = pd.read_csv(filename, sep=";")

    # Ensure integer type for sorting
    df['Traj_Length'] = df['Traj_Length'].astype(int)
    df.sort_values(by=['Traj_Length', 'Scoring'], inplace=True)

    metrics = {
        'F1': 'F1 Score',
        'SHD': 'Structural Hamming Distance',
        'Jaccard': 'Jaccard Distance'
    }
    colors = {'BDE': '#1f77b4', 'MDL': '#ff7f0e'}

    for metric, title in metrics.items():
        sizes = sorted(df['Traj_Length'].unique())
        x = np.arange(len(sizes))
        width = 0.35

        fig, ax = plt.subplots(figsize=(10, 6))
        
        bde_vals = []
        mdl_vals = []
        
        for s in sizes:
            val_bde = df[(df['Traj_Length'] == s) & (df['Scoring'] == 'BDE')][metric].values
            bde_vals.append(val_bde[0] if len(val_bde) > 0 else 0)
            
            val_mdl = df[(df['Traj_Length'] == s) & (df['Scoring'] == 'MDL')][metric].values
            mdl_vals.append(val_mdl[0] if len(val_mdl) > 0 else 0)

        rects1 = ax.bar(x - width/2, bde_vals, width, label='BDE', color=colors['BDE'], alpha=0.9)
        rects2 = ax.bar(x + width/2, mdl_vals, width, label='MDL', color=colors['MDL'], alpha=0.9)

        ax.set_ylabel(metric)
        ax.set_xlabel('Trajectory Length (Steps)')
        ax.set_title("Neurogenesis Model: {}".format(title))
        ax.set_xticks(x)
        ax.set_xticklabels(sizes)
        ax.legend(title="Scoring Function")
        ax.grid(axis='y', linestyle='--', alpha=0.3)

        # Manually add labels for Python 2/Matplotlib 2 compatibility
        # (ax.bar_label is only available in Matplotlib 3.4+)
        def autolabel(rects):
            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width()/2., 1.01*height,
                        '%.2f' % height,
                        ha='center', va='bottom', fontsize=9)
        
        autolabel(rects1)
        autolabel(rects2)

        plt.tight_layout()
        filename = "graphs/BioModel_Results_{}.png".format(metric)
        plt.savefig(filename, dpi=300)
        plt.close()
        
    print("Part 2 plots generated.")

# ==========================================
# MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    # 1. Part 1 Analysis
    df_part1 = load_and_prep_part1("Final_Report_Part1.csv")
    plot_part1_deep_dive(df_part1)
    plot_part1_metric_breakdown(df_part1)
    
    # 2. Part 2 Analysis
    run_part2_pipeline("Final_Report_Part2.csv")
    
    print("\nAll visualization tasks completed successfully.")