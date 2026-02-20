import json
import pandas as pd
import os
import sys
import re

# ==========================================
# CONFIGURATION
# ==========================================
MODELS_FILE = 'results/best_models.json'
METADATA_FILE = 'metadata/metadata.json'
DATA_DIR = 'hypothesis_data'
OUTPUT_REPORT = 'final_model_sample_flags_v6.csv'

# Thresholds for GLOBAL data only
CRITICAL_THRESHOLD_GLOBAL = 5    
WARNING_THRESHOLD_GLOBAL = 30
PROPORTION_THRESHOLD = 0.01   

# Levels that MUST exist in every subgroup (for Likert)
REQUIRED_LIKERT_LEVELS = {1, 2, 3, 4}

def get_group_from_hyp(hyp_name):
    match = re.match(r"^(.*)_\d+$", hyp_name)
    if match:
        return match.group(1)
    return None

def check_global_counts(series, threshold_crit, threshold_warn):
    """
    Checks for LOW SAMPLE SIZE in the entire dataset.
    """
    if series.empty:
        return 'EMPTY', 0, "No data"
    
    min_count = series.min()
    total_n = series.sum()
    proportion = min_count / total_n if total_n > 0 else 0.0
    
    try:
        # For multi-index (interactions), format the tuple nicely
        min_cat = series.idxmin()
    except:
        min_cat = "N/A"

    status = 'OK'
    if min_count < threshold_crit:
        status = 'CRITICAL (Global Low N)'
    elif min_count < threshold_warn:
        status = 'WARNING (Global Low N)'
    
    if status == 'OK' and proportion < PROPORTION_THRESHOLD:
        status = 'WARNING (Global Low %)'

    return status, min_count, f"Smallest level: {min_cat} (n={min_count}, {round(proportion*100,1)}%)"

def check_subgroup_emptiness(df, group_col, var_col, required_levels):
    """
    Checks if any subgroup completely LACKS a required level.
    """
    issues = []
    
    # Get unique values present in each group
    grouped = df.groupby(group_col)[var_col].unique().apply(set)
    
    for group_name, present_levels in grouped.items():
        missing_levels = required_levels - present_levels
        
        if missing_levels:
            issues.append({
                'subgroup': group_name,
                'missing': sorted(list(missing_levels))
            })
            
    return issues

def main():
    print("=== STARTING DIAGNOSTIC CHECK (STRICT: Subgroup Empty + Global Low N + Interactions) ===")
    
    if not os.path.exists(MODELS_FILE) or not os.path.exists(METADATA_FILE):
        print(f"❌ FATAL: Config files missing.")
        return

    with open(MODELS_FILE, 'r') as f: best_models = json.load(f)
    with open(METADATA_FILE, 'r') as f: metadata = json.load(f)

    results = []
    
    # Columns that define subgroups
    GROUPING_COLS = ['topic_condition', 'condition', 'experiment_group', 'topic']

    print("\n--- Processing Models ---")
    for i, (hyp_name, model_info) in enumerate(best_models.items(), 1):
        print(f"[{i}/{len(best_models)}] Checking {hyp_name}...", end=" ")
        
        group_name = get_group_from_hyp(hyp_name)
        if not group_name or group_name not in metadata:
            print(f"❌ SKIPPED (No metadata)")
            continue

        hyp_meta = metadata[group_name].get(hyp_name)
        if not hyp_meta: continue

        dep_var = hyp_meta.get('dep_var')
        
        csv_path = os.path.join(DATA_DIR, group_name, f"{hyp_name}.csv")
        if not os.path.exists(csv_path):
            print(f"❌ CSV Missing")
            continue
            
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"❌ CSV Error")
            continue

        # Identify which grouping column exists in this file
        active_group_col = next((c for c in GROUPING_COLS if c in df.columns), None)

        var_sources = {}
        if dep_var: var_sources[dep_var] = {"Dependent Variable"}
        if 'kept_vars' in model_info:
            for v in model_info['kept_vars']:
                if v not in var_sources: var_sources[v] = set()
                var_sources[v].add("Main Effect")
        
        # -------------------------------------------------
        # 1. CHECK SINGLE VARIABLES
        # -------------------------------------------------
        for var, sources in var_sources.items():
            if var not in df.columns: continue
            
            is_numeric_likert = pd.api.types.is_numeric_dtype(df[var])

            # A. SUBGROUP CHECK: Are there any EMPTY counts? (Strictly for Likert 1-4)
            if active_group_col and is_numeric_likert:
                missing_issues = check_subgroup_emptiness(df, active_group_col, var, REQUIRED_LIKERT_LEVELS)
                for issue in missing_issues:
                    results.append({
                        'Group': group_name, 
                        'Hypothesis': hyp_name, 
                        'Term': var, 
                        'Type': 'Variable', 
                        'Check': f'Emptiness ({active_group_col})', 
                        'Status': 'CRITICAL (Subgroup Empty)', 
                        'Value': 0, 
                        'Details': f"In '{issue['subgroup']}', nobody selected {issue['missing']}."
                    })

            # B. GLOBAL CHECK: Are there Low N counts?
            counts = df[var].value_counts()
            g_status, g_val, g_det = check_global_counts(counts, CRITICAL_THRESHOLD_GLOBAL, WARNING_THRESHOLD_GLOBAL)
            if g_status != 'OK':
                results.append({
                    'Group': group_name, 
                    'Hypothesis': hyp_name, 
                    'Term': var, 
                    'Type': 'Variable', 
                    'Check': 'Sample Size (Global)', 
                    'Status': g_status, 
                    'Value': g_val, 
                    'Details': g_det
                })

        # -------------------------------------------------
        # 2. CHECK INTERACTIONS (Global Only)
        # -------------------------------------------------
        if 'kept_intrs_list' in model_info:
            for intr_vars in model_info['kept_intrs_list']:
                term_name = ":".join(intr_vars)
                
                # Ensure all variables exist
                if not all(v in df.columns for v in intr_vars): continue
                
                # Calculate counts for every combination of the interaction variables
                # groupby().size() gives the count of each unique combination
                intr_counts = df.groupby(intr_vars).size()
                
                # Check Global Counts
                g_status, g_val, g_det = check_global_counts(intr_counts, CRITICAL_THRESHOLD_GLOBAL, WARNING_THRESHOLD_GLOBAL)
                
                if g_status != 'OK':
                    results.append({
                        'Group': group_name, 
                        'Hypothesis': hyp_name, 
                        'Term': term_name, 
                        'Type': 'Interaction', 
                        'Check': 'Sample Size (Global)', 
                        'Status': g_status, 
                        'Value': g_val, 
                        'Details': g_det
                    })

        print("Done.")

    # 3. Save Report
    if results:
        df_res = pd.DataFrame(results)
        
        # Priority Ranking
        status_rank = {
            'CRITICAL (Subgroup Empty)': 0, 
            'CRITICAL (Global Low N)': 1, 
            'WARNING (Global Low N)': 2,
            'WARNING (Global Low %)': 3
        }
        df_res['rank'] = df_res['Status'].map(lambda x: status_rank.get(x, 99))
        df_res = df_res.sort_values(by=['rank', 'Group', 'Hypothesis'])
        df_res = df_res.drop(columns=['rank'])
        
        cols = ['Group', 'Hypothesis', 'Status', 'Check', 'Type', 'Term', 'Value', 'Details']
        df_res = df_res[[c for c in cols if c in df_res.columns]]
        
        df_res.to_csv(OUTPUT_REPORT, index=False)
        print("\n" + "="*60)
        print(f"✅ CHECK COMPLETE. Report saved to: {OUTPUT_REPORT}")
        print("\n--- Top Critical Issues ---")
        print(df_res[df_res['Status'].str.contains('CRITICAL')].head(15).to_string(index=False))
    else:
        print("\n✅ CHECK COMPLETE. No issues found.")

if __name__ == "__main__":
    main()