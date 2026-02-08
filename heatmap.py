import argparse
import itertools
import json
import re
import textwrap
from pathlib import Path
from typing import List, Tuple, Dict, Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib import cm

# ==========================================
# HARDCODED PATHS
# ==========================================
METADATA_PATH = "metadata/metadata.json"  # <--- HARDCODED
CSV_ROOT_DIR  = "hypothesis_data"         # <--- HARDCODED
OUTPUT_DIR    = "heat_map_outputs"                 # <--- HARDCODED

# ==========================================
# CONFIGURATION
# ==========================================
EXCLUDED_COLS = {
    "participant_number", "topic_condition", "usefreq_hf", "usefreq_med",
    "experience_hf", "experience_med", "gender", "ethnicity_combined",
    "age_bracket", "region_broad", "data_practice"
}

PRGn_like = plt.cm.Purples

def sanitize_filename(name: str) -> str:
    return re.sub(r'[\\/*?:"<>|]', "_", str(name))

def _wrap(labels, width=18):
    return ["\n".join(textwrap.wrap(str(x), width=width, break_long_words=False)) for x in labels]

# ==========================================
# PLOTTING FUNCTIONS
# ==========================================
def save_colorbar(fname: Path, vmin: float, vmax: float, label: str):
    fig, ax = plt.subplots(figsize=(1.5, 6))
    norm = Normalize(vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=PRGn_like), cax=ax)
    cbar.set_label(label, fontsize=14, weight='bold')
    cbar.ax.tick_params(labelsize=12)
    fig.subplots_adjust(left=0.1, right=0.3, top=0.9, bottom=0.1)
    fname.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fname.with_suffix(".png"), bbox_inches='tight', dpi=300)
    plt.close(fig)

def plot_heatmap_no_cbar(piv: pd.DataFrame, title: str, fname: Path, fmt: str, vmin: float, vmax: float):
    if piv.size == 0 or piv.notna().sum().sum() == 0: return False

    xt = _wrap(list(piv.columns))
    yt = _wrap(list(piv.index))
    w = max(8, 1.2 * piv.shape[1] + 2)
    h = max(6, 1.0 * piv.shape[0] + 2)

    fig = plt.figure(figsize=(w, h), constrained_layout=True)
    ax = fig.add_subplot(111)
    norm = Normalize(vmin=vmin, vmax=vmax)
    ax.imshow(piv.values, aspect="auto", cmap=PRGn_like, norm=norm, origin="lower")

    ax.set_xticks(np.arange(piv.shape[1])); ax.set_xticklabels(xt, rotation=45, ha="right", fontsize=14)
    ax.set_yticks(np.arange(piv.shape[0])); ax.set_yticklabels(yt, fontsize=14)
    ax.set_xlabel(piv.columns.name, fontsize=18, weight='bold', labelpad=10)
    ax.set_ylabel(piv.index.name, fontsize=18, weight='bold', labelpad=10)

    total_cells = piv.size
    ann_font = 30 if total_cells <= 50 else 16
    if total_cells <= 400:
        midpoint = vmin + (vmax - vmin) * 0.5
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                v = piv.values[i, j]
                if pd.notna(v):
                    txt = fmt.format(v)
                    color = "white" if v > midpoint else "black"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=ann_font, color=color, weight='bold')

    fname.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fname.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True

# ==========================================
# DATA PROCESSING
# ==========================================
def process_frequencies(df: pd.DataFrame, ind_vars: List[str], output_dir: Path):
    valid_cols = [c for c in ind_vars if c in df.columns and c not in EXCLUDED_COLS]
    valid_cols = list(dict.fromkeys(valid_cols)) # remove dupes
    if len(valid_cols) < 2: return

    pairs = list(itertools.combinations(valid_cols, 2))
    freq_dir = output_dir / "frequencies"
    freq_dir.mkdir(parents=True, exist_ok=True)
    
    plot_data = []
    global_max = 0.0

    for col_a, col_b in pairs:
        ct = pd.crosstab(df[col_b], df[col_a], dropna=True)
        if ct.size < 2: continue
        
        total_n = ct.sum().sum()
        if total_n == 0: continue
        
        ct_pct = (ct / total_n) * 100
        global_max = max(global_max, ct_pct.max().max())
        
        plot_data.append({
            "piv": ct_pct,
            "title": f"Frequency: {col_b} vs {col_a}",
            "fname": freq_dir / f"{sanitize_filename(col_b)}_vs_{sanitize_filename(col_a)}"
        })
    
    if not plot_data: return
    if global_max == 0: global_max = 100.0

    save_colorbar(freq_dir / "_Scale_Legend", 0, global_max, "Percentage (%)")
    for item in plot_data:
        plot_heatmap_no_cbar(item["piv"], item["title"], item["fname"], "{:.1f}%", 0, global_max)
    print(f"      -> Generated {len(plot_data)} frequency maps.")

def process_averages(df: pd.DataFrame, dep_var: str, interactions: List[Any], output_dir: Path):
    if dep_var not in df.columns: return

    valid_pairs = []
    for group in interactions:
        if isinstance(group, list) and len(group) == 2:
            if group[0] in df.columns and group[1] in df.columns:
                valid_pairs.append((group[0], group[1]))
    if not valid_pairs: return

    try: numeric_target = pd.to_numeric(df[dep_var], errors='coerce')
    except: return
    if numeric_target.notna().sum() == 0: return

    avg_dir = output_dir / f"avg_{sanitize_filename(dep_var)}"
    avg_dir.mkdir(parents=True, exist_ok=True)

    obs_min, obs_max = numeric_target.min(), numeric_target.max()
    if obs_min == obs_max: obs_min -= 0.5; obs_max += 0.5

    plot_data = []
    for col_a, col_b in valid_pairs:
        temp = df.copy()
        temp['target'] = numeric_target
        grouped = temp.groupby([col_b, col_a])['target'].mean()
        if grouped.empty: continue
        
        piv = grouped.unstack(level=1)
        piv.index.name, piv.columns.name = col_b, col_a
        plot_data.append({
            "piv": piv,
            "title": f"Avg {dep_var}: {col_b} vs {col_a}",
            "fname": avg_dir / f"{sanitize_filename(col_b)}_vs_{sanitize_filename(col_a)}"
        })

    if not plot_data: return
    save_colorbar(avg_dir / "_Scale_Legend", obs_min, obs_max, f"Avg {dep_var}")
    for item in plot_data:
        plot_heatmap_no_cbar(item["piv"], item["title"], item["fname"], "{:.2f}", obs_min, obs_max)
    print(f"      -> Generated {len(plot_data)} average maps for {dep_var}.")

# ==========================================
# MAIN EXECUTION
# ==========================================
def main():
    # Setup Paths
    meta_path = Path(METADATA_PATH)
    csv_root = Path(CSV_ROOT_DIR)
    out_root = Path(OUTPUT_DIR)

    if not meta_path.exists():
        print(f"Error: Metadata file not found at {meta_path}")
        return

    with open(meta_path, 'r') as f:
        raw_metadata = json.load(f)

    # Recursive walker to handle nested structure (e.g. v1 -> v1_0)
    def process_node(node, current_path_parts):
        for key, value in node.items():
            if isinstance(value, dict):
                # Check if this is a leaf configuration node
                if "ind_vars" in value or "dep_var" in value:
                    hypothesis_id = key
                    
                    # Construct CSV path based on nesting
                    # e.g. hypothesis_data/v1/v1_0.csv
                    subfolder = Path(*current_path_parts)
                    csv_file = csv_root / subfolder / f"{hypothesis_id}.csv"
                    
                    print(f"\n--- Processing: {hypothesis_id} ---")
                    
                    if csv_file.exists():
                        try:
                            df = pd.read_csv(csv_file)
                            
                            # Output dir matches structure
                            hyp_out_dir = out_root / subfolder / hypothesis_id
                            hyp_out_dir.mkdir(parents=True, exist_ok=True)
                            
                            # Frequencies
                            if "ind_vars" in value:
                                process_frequencies(df, value["ind_vars"], hyp_out_dir)
                            
                            # Averages
                            if "dep_var" in value and "interactions" in value:
                                process_averages(df, value["dep_var"], value["interactions"], hyp_out_dir)
                                
                        except Exception as e:
                            print(f"    [Error] Failed to process CSV: {e}")
                    else:
                        print(f"    [Skipping] File not found: {csv_file}")
                else:
                    # It's a folder/grouping node, recurse deeper
                    new_path = current_path_parts + [key]
                    process_node(value, new_path)

    print(f"Scanning metadata from {meta_path}...")
    print(f"Looking for CSVs in {csv_root}...")
    process_node(raw_metadata, [])
    print(f"\nDone. Results saved to: {out_root.resolve()}")

if __name__ == "__main__":
    main()