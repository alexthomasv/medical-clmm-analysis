import os
import glob
import pandas as pd
import json
import re

# ==============================================================================
# CONFIGURATION
# ==============================================================================
RESULTS_DIR = "results"
OUTPUT_JSON = "best_models.json"

# ==============================================================================
# HELPERS
# ==============================================================================

def get_all_run_dirs(base_path):
    if not os.path.exists(base_path):
        return []
    return sorted([d for d in glob.glob(os.path.join(base_path, "run_*")) if os.path.isdir(d)])

def parse_interaction_string(interaction_str):
    """
    Parses 'scale(as.numeric(x)):scale(as.numeric(y))' or 'x:y'
    Returns a list of parent variables involved.
    """
    if not interaction_str or pd.isna(interaction_str):
        return None
    
    # Split by colon to get parts
    parts = interaction_str.split(":")
    
    # Clean up R syntax wrappers
    clean_parts = []
    for p in parts:
        # Remove scale(as.numeric(...)) wrapper
        p = re.sub(r"scale\(as\.numeric\(`?", "", p)
        p = re.sub(r"`?\)\)", "", p)
        # Remove just backticks if present
        p = p.replace("`", "").strip()
        clean_parts.append(p)
        
    return clean_parts

# ==============================================================================
# MAIN LOGIC
# ==============================================================================

def validate_and_generate():
    print(f"Scanning '{RESULTS_DIR}' for runs...")
    run_dirs = get_all_run_dirs(RESULTS_DIR)
    
    if not run_dirs:
        print("No runs found.")
        return

    print(f"Found {len(run_dirs)} runs.")
    
    # Data structure: { "hyp_id": { "aic": val, "predictors": "str", "interactions": "str" } }
    consensus_data = {}
    inconsistencies = []

    # 1. READ AND VALIDATE
    for rdir in run_dirs:
        run_name = os.path.basename(rdir)
        csv_path = os.path.join(rdir, "aic_results.csv")
        
        if not os.path.exists(csv_path):
            print(f"Skipping {run_name} (no CSV)")
            continue
            
        try:
            df = pd.read_csv(csv_path)
            df.columns = df.columns.str.strip()
            
            required = ["hypothesis", "AIC", "kept_predictors", "interactions"]
            if not all(col in df.columns for col in required):
                print(f"Skipping {run_name} (Missing columns)")
                continue

            for _, row in df.iterrows():
                hyp = str(row['hypothesis']).strip()
                aic = float(row['AIC']) if pd.notna(row['AIC']) else None
                preds = str(row['kept_predictors']).strip() if pd.notna(row['kept_predictors']) else ""
                intrs = str(row['interactions']).strip() if pd.notna(row['interactions']) else ""

                if aic is None: continue # Skip failed models

                current_data = (aic, preds, intrs)

                if hyp not in consensus_data:
                    consensus_data[hyp] = current_data
                else:
                    # VALIDATION CHECK
                    # We check strict AIC equality (or close tolerance)
                    existing_aic = consensus_data[hyp][0]
                    
                    if abs(existing_aic - aic) > 1e-5:
                        inconsistencies.append((hyp, run_name, existing_aic, aic))

        except Exception as e:
            print(f"Error reading {run_name}: {e}")

    # 2. HANDLE ERRORS
    if inconsistencies:
        print("\n" + "!"*60)
        print(f"FATAL ERROR: AIC Mismatch found in {len(inconsistencies)} hypotheses")
        print("!"*60)
        for hyp, run, old_aic, new_aic in inconsistencies:
            print(f"  {hyp}: Existing={old_aic:.4f} vs Run({run})={new_aic:.4f}")
        
        # Stop execution if validation fails
        assert False, "Validation Failed: Inconsistent results across runs."

    print("\nValidation PASSED. All runs are consistent.")

    # 3. GENERATE JSON FOR R
    print("Generating best_models.json...")
    
    json_output = {}
    
    for hyp, (aic, preds_str, intrs_str) in consensus_data.items():
        
        # Parse Predictors list
        kept_vars = [p.strip() for p in preds_str.split("+") if p.strip()]
        
        # Parse Interactions list
        kept_intrs = [i.strip() for i in intrs_str.split("+") if i.strip()]
        
        # Reconstruct Parents structure (Required by R script)
        interaction_parents = []
        for i_str in kept_intrs:
            parents = parse_interaction_string(i_str)
            if parents:
                interaction_parents.append(parents)

        json_output[hyp] = {
            "kept_vars": kept_vars,
            "kept_intrs": kept_intrs,
            "interaction_parents": interaction_parents
        }

    # Write file
    out_path = os.path.join(RESULTS_DIR, OUTPUT_JSON)
    with open(out_path, 'w') as f:
        json.dump(json_output, f, indent=4)

    print(f"Success! JSON saved to: {out_path}")
    print(f"Contains {len(json_output)} models ready for final fitting.")

if __name__ == "__main__":
    validate_and_generate()