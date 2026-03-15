# Loads v1/v2 hypothesis

from scipy.stats import kstest, norm, kruskal 
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant

import gspread
import google.auth

import re
from collections import defaultdict
import pandas as pd
import numpy as np 

from logging import exception
from functools import reduce
import os
import json
from pathlib import Path
from scipy.stats import levene

# ==========================================
# 1. UTILITY FUNCTIONS
# ==========================================

def get_nickname(string):
  pattern = r"_\d$"
  extracted_string = re.sub(pattern, "", string)
  return extracted_string

def parse_float_maybe(string_val):
    cleaned = re.sub(r"\(.*?\)", "", string_val)
    cleaned = re.sub(r"[^\d.+\-eE]", "", cleaned)
    if not cleaned:
        return None
    if not re.match(r"^[+\-]?\d+(\.\d+)?([eE][+\-]?\d+)?$", cleaned):
        return None
    return float(cleaned)

def flatten_summary_dict(nested_dict):
  flat_dict = {}
  def flatten(current_key, value):
      if isinstance(value, dict):
          for k, v in value.items():
              new_key = f"{current_key} - {k}" if current_key else k
              flatten(new_key, v)
      else:
          flat_dict[current_key] = value
  for key, val in nested_dict.items():
      flatten(key, val)
  return flat_dict

def has_paren(s: str) -> bool:
    return re.search(r'\(\s*.*?\s*\)', s) is not None

def recode_experience(val):
    s = str(val).strip()
    if s in ['1', '2']:
        return 'Baddish'
    elif s in ['3', '4']:
        return 'Goodish'
    elif 'never' in s.lower():
        return 'Never Used One'
    return s 

def recode_ethnicity_combined(val: str) -> str:
  s = str(val).strip()
  if s in ['NO_DATA']:
      return 'Other'
  return s

# Helper to parse a formula string like "var1 + var2 + var1*var2"
def parse_formula_string(formula_str):
    if not formula_str:
        return [], []
    
    parts = formula_str.split("+")
    main_vars = []
    interactions = []
    
    for x in parts:
        x = x.strip()
        if not x or has_paren(x): # Skip empty or random effects (1|...)
            continue
            
        if '*' in x:
            # It's an interaction (e.g. A*B)
            terms = [t.strip() for t in x.split('*')]
            interactions.append(tuple(terms))
            # Ensure constituents are in main vars list for data pulling
            for t in terms:
                if t not in main_vars:
                    main_vars.append(t)
        elif ':' in x:
             # R style interaction (A:B)
            terms = [t.strip() for t in x.split(':')]
            interactions.append(tuple(terms))
            for t in terms:
                if t not in main_vars:
                    main_vars.append(t)
        else:
            if x not in main_vars:
                main_vars.append(x)
                
    return main_vars, interactions

# ==========================================
# 2. HYPOTHESIS CLASS
# ==========================================

class Hypothesis:
  def __init__(self, num, dep_var, 
               all_data_cols, # All columns needed to pull data
               base_preds, base_inters, # The Baseline Model
               cand_preds, cand_inters, # The Candidates (v6 only)
               app_topics):
    
    self.num = num
    self.dep_vars = [dep_var]
    self.ind_vars = all_data_cols # Used for creating the CSV
    
    # Metadata Fields
    self.base_predictors = base_preds
    self.base_interactions = base_inters
    self.candidate_predictors = cand_preds
    self.candidate_interactions = cand_inters
    
    self.app_topics = app_topics

  def get_interaction_strings(self, inter_list):
      """Converts list of tuples [('A','B')] to list of strings ['A:B']"""
      return [":".join(i) for i in inter_list]

  def get_base_predictors_combined(self):
      """Helper for old code logic if needed"""
      filter_set = {'usefreq_hf', 'experience_hf', 'gender', 'ethnicity_combined', 'age_bracket', 'region_broad'}
      base_set = []
      for v in self.base_predictors:
          if v not in filter_set: base_set.append(v)
      for inter in self.base_interactions:
          base_set.append(":".join(inter))
      return list(set(base_set))

  def __repr__(self):
    return f"{self.dep_vars[0]} ({self.num})"

# ==========================================
# 3. DATA EXTRACTION
# ==========================================

def extract_hypothesis_special(hypothesis_sheet, hypothesis_sheet_name, app_topics):
  hypothesis_list = []
  all_rows = hypothesis_sheet.get_all_values()
  headers = [h.lower().strip() for h in all_rows[0]] 
  
  # 1. Identify Columns
  try:
      idx_dep = headers.index("dv")
      idx_baseline = headers.index("baseline model")
  except ValueError as e:
      print(f"Error processing {hypothesis_sheet_name}: Missing required column. {e}")
      return []

  # Check if v6 (Forward Selection Mode)
  is_v6 = "v6" in hypothesis_sheet_name
  idx_all = -1
  if is_v6:
      try:
          idx_all = headers.index("all predictors")
      except ValueError:
          print(f"Error: v6 sheet '{hypothesis_sheet_name}' requires 'All Predictors' column.")
          return []

  h_idx = 0
  for row in all_rows[1:]:
    if not row or not row[0]: continue 

    dep_var = row[idx_dep].strip()
    
    # --- A. Parse Baseline Model (Everyone has this) ---
    baseline_str = row[idx_baseline]
    base_vars, base_inters = parse_formula_string(baseline_str)
    
    # Initialize Candidates
    cand_vars = []
    cand_inters = []
    
    # --- B. Parse Candidates (Only v6) ---
    if is_v6:
        all_preds_str = row[idx_all]
        all_vars, all_inters = parse_formula_string(all_preds_str)
        
        # Calculate Difference: All - Baseline = Candidates
        cand_vars = list(set(all_vars) - set(base_vars))
        cand_inters = list(set(all_inters) - set(base_inters))
    
    # --- C. Consolidate Data Columns ---
    # We need to pull data for EVERYTHING (Base + Candidates)
    all_data_cols = list(set(base_vars + cand_vars))
    
    # Create Object
    h_num = f"{hypothesis_sheet_name}_{h_idx}"
    hyp = Hypothesis(h_num, dep_var, all_data_cols, 
                     base_vars, base_inters,
                     cand_vars, cand_inters, 
                     app_topics)
    hypothesis_list.append(hyp)
    h_idx += 1

  return hypothesis_list

def extract_predictors(data):
  predictor_map = {}
  data_rows = data.get_all_values()
  predictor_to_data = defaultdict(set)
  predictor_names = data_rows[0][2:]
  data_practices = data_rows[1][2:]
  
  for data_row in data_rows[2:]:
    if len(data_row) < 1509: 
        print(f"Skipping row {data_row} because it has less than 1509 columns")
        continue

    participant_number = data_row[0]
    app_topic = data_row[1]
    usefreq_hf           = data_row[1494]
    usefreq_med          = data_row[1495]
    experience_hf        = recode_experience(data_row[1496])
    experience_med       = recode_experience(data_row[1497])
    gender               = data_row[1498]
    ethnicity_combined   = recode_ethnicity_combined(data_row[1502])
    age_bracket          = data_row[1504]
    region_broad         = data_row[1507]

    assert len(data_row[2:]) == len(data_practices), f"Data row {data_row} has {len(data_row)} columns, but {len(data_practices)} data practices"
    for predictor, data_practice, val in zip(predictor_names, data_practices, data_row[2:]):
        if not val: 
            continue
            
        # --- FIX: Strip accidental whitespace from the practice string ---
        clean_practice = data_practice.strip()
        
        predictor_to_data[predictor].add((
            participant_number, app_topic, usefreq_hf, usefreq_med, 
            experience_hf, experience_med, gender, ethnicity_combined, 
            age_bracket, region_broad, clean_practice, val
        ))
        
  for predictor in predictor_to_data:
    data = list(predictor_to_data[predictor])
    cols = ['participant_number', 'topic_condition', 'usefreq_hf', 'usefreq_med', 'experience_hf', 'experience_med', 'gender', 'ethnicity_combined', 'age_bracket', 'region_broad', 'data_practice', predictor]
    if cols.count(predictor) > 1: cols[-1] = f"{predictor}_target"

    df = pd.DataFrame(data, columns=cols)
    valid_cats = [x for x in ["Goodish", "Baddish", "Never Used One"] if x in df['experience_hf'].unique()]
    df['experience_hf'] = pd.Categorical(df['experience_hf'], categories=valid_cats, ordered=False)
    df_sorted = df.sort_values(by=['data_practice',  'topic_condition', 'participant_number'])
    predictor_map[predictor] = df_sorted

  return predictor_map

class Survey:
  def __init__(self, predictor_map):
    self.predictor_map = predictor_map
    self.missing_data_log = [] # NEW: Initialize a master list to track missing data

  def create_local_table(self, hypothesis: Hypothesis, topic_conditions_filter):
    exceptions = ['purpose_health', 'purpose_sens']
    special = ['region_broad', 'usefreq_med', 'experience_hf', 'experience_med', 'gender', 'ethnicity_combined', 'usefreq_hf', 'age_bracket']
    
    # Identify all needed columns from ind_vars (which includes base + candidates)
    needed_vars = set(hypothesis.base_predictors + hypothesis.candidate_predictors)
    for inter in hypothesis.base_interactions + hypothesis.candidate_interactions:
        for component in inter:
            needed_vars.add(component)
    
    ind_vars = list(needed_vars)
    dep_var = hypothesis.dep_vars[0]

    # Filter out special cols from predictor list
    exception_predictors = []
    filtered_predictors = []

    for var in ind_vars:
      if var in special:
        continue
      if var in exceptions:
        exception_predictors.append(var)
      else:
        filtered_predictors.append(var)

    ind_vars_df = [self.predictor_map[ind_var] for ind_var in filtered_predictors]
    exception_predictors_df = [self.predictor_map[exception].drop('data_practice', axis=1) for exception in exception_predictors]

    dep_var_df = self.predictor_map[dep_var]
    
    # --- DETAILED MISSING DATA REPORTING LOGIC ---
    columns = ind_vars_df + [dep_var_df]
    var_names = filtered_predictors + [dep_var] # Track the actual variable names
    
    merge_keys = ['participant_number', 'topic_condition', 'data_practice']
    sets_of_keys = []
    
    for df in columns:
        keys = set(df[merge_keys].itertuples(index=False, name=None))
        sets_of_keys.append(keys)
        
    if sets_of_keys:
        all_keys = set.union(*sets_of_keys)
        surviving_keys = set.intersection(*sets_of_keys)
        dropped_keys = all_keys - surviving_keys
        
        if dropped_keys:
            print(f"\n[!] Missing Data for Hypothesis {hypothesis.num}: Dropped {len(dropped_keys)} rows.")
            
            for i, key in enumerate(sorted(list(dropped_keys))):
                # Find which sets DO NOT contain this key, and grab the corresponding variable name
                missing_vars = [var_names[idx] for idx, k_set in enumerate(sets_of_keys) if key not in k_set]
                missing_str = ", ".join(missing_vars)
                
                # Append to the master audit log
                self.missing_data_log.append({
                    "Hypothesis": hypothesis.num,
                    "Participant": key[0],
                    "Topic": key[1],
                    "Practice": key[2],
                    "Missing_Variables": missing_str
                })
                
                # Console preview (Capped at 5 so the terminal stays clean)
                if i < 5: 
                    print(f"    - Participant: {key[0]} | Topic: {key[1]} | Practice: '{key[2][:35]}...'")
                    print(f"      -> Missing: {missing_str}")
            
            if len(dropped_keys) > 5:
                print(f"    ... and {len(dropped_keys) - 5} more rows logged to audit file.")
    # -------------------------------------------------

    # First merge block
    final = reduce(lambda left, right: pd.merge(left, right, on=['participant_number', 'topic_condition', 'data_practice'] + special, how='inner'), columns)

    # Second merge block (if exceptions exist)
    if exception_predictors_df:
      final = reduce(
          lambda L, R: pd.merge(L, R, on=['participant_number','topic_condition'] + special, how='inner'),
          [final] + exception_predictors_df
      )
    
    if topic_conditions_filter: 
        final = final[final['topic_condition'].isin(topic_conditions_filter)]

    return final

  def prepare_clmm_data(self, hypothesis_dict):
    os.makedirs("hypothesis_data", exist_ok=True)
    os.makedirs("metadata", exist_ok=True)

    metadata = {}
    for hypothesis_name, hypothesis_list in hypothesis_dict.items():
        metadata[hypothesis_name] = {}
        for hypothesis in hypothesis_list:
            local_table = self.create_local_table(hypothesis, hypothesis.app_topics)
            
            # Save Data
            out_dir = f"hypothesis_data/{hypothesis_name}"
            os.makedirs(out_dir, exist_ok=True)
            local_table.to_csv(f"{out_dir}/{hypothesis.num}.csv", index=False)
            
            metadata[hypothesis_name][hypothesis.num] = {
                "dep_var": hypothesis.dep_vars[0],
                "base_predictors": hypothesis.base_predictors,
                "base_interactions": hypothesis.base_interactions,
                "candidate_terms_predictors": hypothesis.candidate_predictors,
                "candidate_terms_interactions": hypothesis.candidate_interactions
            }
            
    with open("metadata/metadata.json", "w") as f:
        json.dump(metadata, f, indent=4)

    # --- NEW: EXPORT THE MISSING DATA AUDIT LOG ---
    if self.missing_data_log:
        audit_df = pd.DataFrame(self.missing_data_log)
        audit_path = "metadata/missing_data_audit.csv"
        audit_df.to_csv(audit_path, index=False)
        print(f"\n[+] Full missing data audit saved to '{audit_path}'")
    else:
        print("\n[+] No missing data detected across any hypothesis!")

# ==========================================
# 4. STATS REPORTER
# ==========================================
class StatsReporter:
    def __init__(self, predictor_map):
        self.predictor_map = predictor_map

    def calculate_stats(self):
        all_records = []
        for predictor_name, df in self.predictor_map.items():
            target_col = f"{predictor_name}_target" if f"{predictor_name}_target" in df.columns else predictor_name
            temp_df = df[['topic_condition', target_col]].copy()
            temp_df.columns = ['App_Topic', 'Value']
            temp_df['Data_Type'] = predictor_name
            temp_df['Value'] = pd.to_numeric(temp_df['Value'], errors='coerce')
            temp_df = temp_df.dropna(subset=['Value'])
            all_records.append(temp_df)
            
        if not all_records:
            return pd.DataFrame()

        master_df = pd.concat(all_records, ignore_index=True)
        report_data = []

        for dtype, group in master_df.groupby('Data_Type'):
            stats = self._compute_stats(group['Value'])
            stats.update({'Level': 'All Apps Per Data Type', 'Data Type': dtype, 'App Topic': 'ALL'})
            report_data.append(stats)

        for (dtype, app), group in master_df.groupby(['Data_Type', 'App_Topic']):
            stats = self._compute_stats(group['Value'])
            stats.update({'Level': 'Data Type Per App', 'Data Type': dtype, 'App Topic': app})
            report_data.append(stats)

        report_df = pd.DataFrame(report_data)
        report_df = report_df.sort_values(by=['Level', 'Data Type', 'App Topic'])
        cols = ['Level', 'Data Type', 'App Topic', 'Mean', 'Median', 'Mode', 'Std Dev', 'Count']
        cols = [c for c in cols if c in report_df.columns]
        return report_df[cols]

    def _compute_stats(self, series):
        if series.empty:
             return {'Mean': np.nan, 'Median': np.nan, 'Mode': np.nan, 'Std Dev': np.nan, 'Count': 0}
        desc = series.describe()
        mode_val = series.mode()
        return {
            'Mean': round(desc['mean'], 2),
            'Median': round(desc['50%'], 2),
            'Mode': mode_val.iloc[0] if not mode_val.empty else np.nan,
            'Std Dev': round(desc['std'], 2),
            'Count': int(desc['count'])
        }

# ==========================================
# 5. STATISTICAL TESTER (KRUSKAL-WALLIS)
# ==========================================
class StatTester:
    def __init__(self, predictor_map):
        self.predictor_map = predictor_map

    def check_app_differences(self, variables_to_test=None):
        results = []
        if variables_to_test is None:
            variables_to_test = list(self.predictor_map.keys())

        print(f"\n--- Running Kruskal-Wallis Tests (Grouping by App Topic) ---")
        for var_name in variables_to_test:
            if var_name not in self.predictor_map:
                continue

            df = self.predictor_map[var_name]
            target_col = f"{var_name}_target" if f"{var_name}_target" in df.columns else var_name
            df = df.copy()
            df[target_col] = pd.to_numeric(df[target_col], errors='coerce')
            df = df.dropna(subset=[target_col])

            groups = [group[target_col].values for name, group in df.groupby('topic_condition')]

            if len(groups) < 2:
                continue

            try:
                stat, p_value = kruskal(*groups)
                significance = ""
                if p_value < 0.001: significance = "***"
                elif p_value < 0.01: significance = "**"
                elif p_value < 0.05: significance = "*"

                results.append({
                    'Variable': var_name,
                    'H-statistic': round(stat, 3),
                    'p-value': p_value, 
                    'Significance': significance
                })
            except Exception as e:
                print(f"Error testing {var_name}: {e}")

        res_df = pd.DataFrame(results)
        if not res_df.empty:
            res_df = res_df.sort_values(by='p-value')
        return res_df

# ==========================================
# 6. VIF ANALYZER
# ==========================================
class VIFAnalyzer:
    def __init__(self, survey_instance):
        self.survey = survey_instance

    def compute_vif_for_hypothesis(self, df, ind_vars):
        valid_predictors = [v for v in ind_vars if v in df.columns]
        if not valid_predictors: return pd.DataFrame()

        X_processed = pd.DataFrame()
        var_groups = {} # Maps original predictor name to its generated column(s)
        
        # 1. Map Ordinal/Likert scales to continuous numeric values
        freq_map = {
            "Less than once a year": 1,
            "A few times a year": 2,
            "A few times a month": 3,
            "A few times a week": 4,
            "Daily": 5
        }
        exp_map = {
            "Never Used One": 1,
            "Baddish": 2,
            "Goodish": 3
        }
        age_map = {
            "18-20": 1,
            "21-44": 2,
            "45-64": 3,
            "65+": 4
        }

        # 2. Process predictors: map ordinals, dummy-encode nominals
        for var in valid_predictors:
            if 'usefreq' in var:
                X_processed[var] = pd.to_numeric(df[var].replace(freq_map), errors='coerce')
                var_groups[var] = [var]
            elif 'experience' in var:
                X_processed[var] = pd.to_numeric(df[var].replace(exp_map), errors='coerce')
                var_groups[var] = [var]
            elif 'age_bracket' in var:
                X_processed[var] = pd.to_numeric(df[var].replace(age_map), errors='coerce')
                var_groups[var] = [var]
            else:
                # Nominal categories (e.g., gender, region) get converted to dummies
                if df[var].dtype == 'object' or isinstance(df[var].dtype, pd.CategoricalDtype):
                    dummies = pd.get_dummies(df[var], prefix=var, drop_first=True, dtype=float)
                    for col in dummies.columns:
                        X_processed[col] = dummies[col]
                    var_groups[var] = list(dummies.columns) # Group them together!
                else:
                    # Fallback for any standard numeric data
                    X_processed[var] = pd.to_numeric(df[var], errors='coerce')
                    var_groups[var] = [var]

        # Drop any rows with missing data after processing
        X_processed = X_processed.dropna()
        if X_processed.empty or X_processed.shape[1] < 2: 
            return pd.DataFrame()

        # 3. Calculate VIF / GVIF
        # We use the correlation matrix determinant method for GVIF
        corr_matrix = X_processed.corr().values
        det_R = np.linalg.det(corr_matrix)

        vif_data = []
        for var_name, cols in var_groups.items():
            if not cols:
                continue
            
            # A. Standard VIF for single-column predictors (Likert scales)
            if len(cols) == 1:
                col_idx = X_processed.columns.get_loc(cols[0])
                X_with_const = sm.add_constant(X_processed.values, has_constant='add')
                try:
                    val = variance_inflation_factor(X_with_const, col_idx + 1)
                    vif_data.append({'Predictor': var_name, 'VIF': round(val, 4)})
                except Exception as e:
                    vif_data.append({'Predictor': var_name, 'VIF': str(e)})
            
            # B. Generalized VIF (GVIF) for multi-column nominal predictors
            else:
                col_indices = [X_processed.columns.get_loc(c) for c in cols]
                comp_indices = [i for i in range(X_processed.shape[1]) if i not in col_indices]
                
                try:
                    # Calculate determinants for the specific variable group vs everything else
                    det_R11 = np.linalg.det(corr_matrix[np.ix_(col_indices, col_indices)])
                    det_R22 = np.linalg.det(corr_matrix[np.ix_(comp_indices, comp_indices)])
                    
                    if det_R == 0:
                        gvif = float('inf') # Perfect multicollinearity
                    else:
                        gvif = (det_R11 * det_R22) / det_R
                        
                    vif_data.append({'Predictor': var_name, 'VIF': round(gvif, 4)})
                except Exception as e:
                    vif_data.append({'Predictor': var_name, 'VIF': str(e)})

        return pd.DataFrame(vif_data)

    def analyze_all_hypotheses(self, hypothesis_dict):
        results = []
        print("\n--- Starting VIF Analysis per Hypothesis ---")
        
        for group_name, h_list in hypothesis_dict.items():
            for h in h_list:
                df = self.survey.create_local_table(h, h.app_topics)
                
                vars_to_check = h.base_predictors 
                if not vars_to_check: vars_to_check = h.ind_vars 
                
                vif_df = self.compute_vif_for_hypothesis(df, vars_to_check)
                
                if not vif_df.empty:
                    vif_df['Hypothesis Group'] = group_name
                    vif_df['Hypothesis Num'] = h.num
                    cols = ['Hypothesis Group', 'Hypothesis Num', 'Predictor', 'VIF']
                    results.append(vif_df[cols])
        
        if results: return pd.concat(results, ignore_index=True)
        return pd.DataFrame()

def perform_levenes_test(predictor_map):
    """
    Performs Levene's test for equal variances across the Low, Medium, 
    and High healthishness app bands, filtered for high-healthishness practices.
    Also includes an unfiltered test for purpose_health.
    """
    print("\n" + "="*40)
    print("PERFORMING LEVENE'S TEST FOR EQUAL VARIANCES")
    print("="*40)
    
    # 1. Define App Groups based on your existing topics
    app_groups = {
        "Low Health Apps": ["Flashlight", "Face Editor", "Kinky Dating"],
        "Medium Health Apps": ["STI Dating", "Marijuana Shopping", "Blood Donation"],
        "High Health Apps": ["Brain Exercises", "Step Counter", "Smart Scale", "STI Help", "Period Tracker"]
    }
    
    # --- ASSERTIONS: Ensure base healthishness metrics exist ---
    assert 'dtype_health' in predictor_map, "Fatal Error: 'dtype_health' not found in predictor_map."
    assert 'duse_health' in predictor_map, "Fatal Error: 'duse_health' not found in predictor_map."
    assert 'recipient_health' in predictor_map, "Fatal Error: 'recipient_health' not found in predictor_map."

    # 2. Identify the High-Healthishness Practices (>= 3.0)
    _, high_dtype_df = find_high_average_practices(predictor_map, 'dtype_health', 3.0)
    high_dtypes = high_dtype_df['data_practice'].tolist() if not high_dtype_df.empty else []
    
    _, high_duse_df = find_high_average_practices(predictor_map, 'duse_health', 3.0)
    high_duses = high_duse_df['data_practice'].tolist() if not high_duse_df.empty else []
    
    _, high_recip_df = find_high_average_practices(predictor_map, 'recipient_health', 3.0)
    high_recipients = high_recip_df['data_practice'].tolist() if not high_recip_df.empty else []
        
    # 3. Configure the metrics to analyze
    metrics_config = [
        ("Dtype Metrics", ['dtype_likely', 'do_strules_dtype', 'should_strules_dtype'], high_dtypes),
        ("Duse Metrics", ['duse_likely', 'do_strules_duse', 'should_strules_duse'], high_duses),
        ("Recipient Metrics", ['recipient_likely', 'do_strules_recipient', 'should_strules_recipient'], high_recipients)
    ]
    
    results = []
    
    # 4. Extract data and perform the test for Dtypes, Duses, and Recipients
    for category_name, metrics, high_practices in metrics_config:
        for metric in metrics:
            if metric not in predictor_map:
                print(f"Skipping {metric}, not found in predictor_map.")
                continue
                
            df = predictor_map[metric].copy()
            target_col = f"{metric}_target" if f"{metric}_target" in df.columns else metric
            
            df[target_col] = pd.to_numeric(df[target_col], errors='coerce')
            df = df.dropna(subset=[target_col, 'topic_condition', 'data_practice'])
            
            # Filter: Only keep High-Health Practices
            df_filtered = df[df['data_practice'].isin(high_practices)]
            
            # Create three separate arrays for the Low, Medium, and High app bands
            low_vals = df_filtered[df_filtered['topic_condition'].isin(app_groups["Low Health Apps"])][target_col].values
            med_vals = df_filtered[df_filtered['topic_condition'].isin(app_groups["Medium Health Apps"])][target_col].values
            high_vals = df_filtered[df_filtered['topic_condition'].isin(app_groups["High Health Apps"])][target_col].values
            
            # Ensure we have data in at least two bands to compare variances
            valid_groups = [g for g in [low_vals, med_vals, high_vals] if len(g) > 1]
            
            if len(valid_groups) >= 2:
                # Run Levene's Test (default center='median' is perfect for Likert scales)
                stat, p_val = levene(*valid_groups, center='median')
                
                sig = ""
                if p_val < 0.001: sig = "***"
                elif p_val < 0.01: sig = "**"
                elif p_val < 0.05: sig = "*"
                
                results.append({
                    "Data_Category": category_name,
                    "Metric": metric,
                    "Levene_Statistic": round(stat, 4),
                    "p-value": p_val,
                    "Significance": sig,
                    "N_Low_Apps": len(low_vals),
                    "N_Med_Apps": len(med_vals),
                    "N_High_Apps": len(high_vals)
                })
            else:
                results.append({
                    "Data_Category": category_name,
                    "Metric": metric,
                    "Levene_Statistic": np.nan,
                    "p-value": np.nan,
                    "Significance": "Not enough data",
                    "N_Low_Apps": len(low_vals),
                    "N_Med_Apps": len(med_vals),
                    "N_High_Apps": len(high_vals)
                })
                
    # 5. SPECIAL CASE: Unfiltered Levene's test for purpose_health
    assert 'purpose_health' in predictor_map, "Fatal Error: 'purpose_health' not found in predictor_map."
    df_purpose = predictor_map['purpose_health'].copy()
    target_col_purpose = "purpose_health_target" if "purpose_health_target" in df_purpose.columns else "purpose_health"
    
    df_purpose[target_col_purpose] = pd.to_numeric(df_purpose[target_col_purpose], errors='coerce')
    # We only drop NAs on the target and topic condition, completely ignoring data_practice constraints
    df_purpose = df_purpose.dropna(subset=[target_col_purpose, 'topic_condition'])
    
    low_vals_p = df_purpose[df_purpose['topic_condition'].isin(app_groups["Low Health Apps"])][target_col_purpose].values
    med_vals_p = df_purpose[df_purpose['topic_condition'].isin(app_groups["Medium Health Apps"])][target_col_purpose].values
    high_vals_p = df_purpose[df_purpose['topic_condition'].isin(app_groups["High Health Apps"])][target_col_purpose].values
    
    valid_groups_p = [g for g in [low_vals_p, med_vals_p, high_vals_p] if len(g) > 1]
    
    if len(valid_groups_p) >= 2:
        stat_p, p_val_p = levene(*valid_groups_p, center='median')
        
        sig_p = ""
        if p_val_p < 0.001: sig_p = "***"
        elif p_val_p < 0.01: sig_p = "**"
        elif p_val_p < 0.05: sig_p = "*"
        
        results.append({
            "Data_Category": "Purpose Metrics",
            "Metric": "purpose_health",
            "Levene_Statistic": round(stat_p, 4),
            "p-value": p_val_p,
            "Significance": sig_p,
            "N_Low_Apps": len(low_vals_p),
            "N_Med_Apps": len(med_vals_p),
            "N_High_Apps": len(high_vals_p)
        })
    else:
        assert False, "No valid groups found for purpose_health."
                
    return pd.DataFrame(results)

# ==========================================
# 7. METRIC AND POSTERITY FUNCTIONS
# ==========================================

def find_high_average_practices(predictor_map, variable_name='dtype_health', threshold=3.0):
    """
    Calculates the average, median, mode, stdev, and variance of a specific variable 
    for each data practice, and returns both the full stats list and a filtered list 
    of practices meeting the average threshold.
    """
    if variable_name not in predictor_map:
        assert False, f"{variable_name} not found in predictor_map."

    # Get the dataframe for this specific variable
    df = predictor_map[variable_name].copy()
    
    # Identify the target column name
    target_col = f"{variable_name}_target" if f"{variable_name}_target" in df.columns else variable_name
    
    # Convert to numeric just in case
    df[target_col] = pd.to_numeric(df[target_col], errors='coerce')
    df = df.dropna(subset=[target_col, 'data_practice'])

    # Helper function to extract mode safely
    def get_mode(x):
        m = x.mode()
        return m.iloc[0] if not m.empty else np.nan

    # 1. Calculate ALL stats at once
    stats_df = df.groupby('data_practice')[target_col].agg(
        Average='mean',
        Median='median',
        Mode=get_mode,
        Stdev='std',
        Variance='var',
        Count='count'
    ).reset_index()

    # Rename 'Average' so it matches your original output format
    stats_df.rename(columns={'Average': f'Average_{variable_name}'}, inplace=True)

    # Round float columns for cleaner CSV output
    float_cols = [f'Average_{variable_name}', 'Stdev', 'Variance']
    stats_df[float_cols] = stats_df[float_cols].round(4)

    # 2. Filter for those with an average >= threshold
    high_avg_df = stats_df[stats_df[f'Average_{variable_name}'] >= threshold].copy()
    
    # Sort them from highest to lowest for easy reading
    stats_df = stats_df.sort_values(by=f'Average_{variable_name}', ascending=False)
    high_avg_df = high_avg_df.sort_values(by=f'Average_{variable_name}', ascending=False)

    return stats_df, high_avg_df

def calculate_stdevs(predictor_map):
    """
    Groups apps by Low/Medium/High healthishness and calculates the standard 
    deviations for specific metrics, filtering only for practices (dtypes, 
    duses, recipients) that have an average healthishness >= 3.0.
    """
    # 1. Define App Groups based on your existing topics
    app_groups = {
        "Low Health Apps": ["Flashlight", "Face Editor", "Kinky Dating"],
        "Medium Health Apps": ["STI Dating", "Marijuana Shopping", "Blood Donation"],
        "High Health Apps": ["Brain Exercises", "Step Counter", "Smart Scale", "STI Help", "Period Tracker"]
    }
    
    # --- ASSERTIONS: Ensure base healthishness metrics exist ---
    assert 'dtype_health' in predictor_map, "Fatal Error: 'dtype_health' not found in predictor_map."
    assert 'duse_health' in predictor_map, "Fatal Error: 'duse_health' not found in predictor_map."
    assert 'recipient_health' in predictor_map, "Fatal Error: 'recipient_health' not found in predictor_map."

    # 2. Identify the High-Healthishness Practices (>= 3.0)
    _, high_dtype_df = find_high_average_practices(predictor_map, 'dtype_health', 3.0)
    high_dtypes = high_dtype_df['data_practice'].tolist() if not high_dtype_df.empty else []
    
    _, high_duse_df = find_high_average_practices(predictor_map, 'duse_health', 3.0)
    high_duses = high_duse_df['data_practice'].tolist() if not high_duse_df.empty else []
    
    _, high_recip_df = find_high_average_practices(predictor_map, 'recipient_health', 3.0)
    high_recipients = high_recip_df['data_practice'].tolist() if not high_recip_df.empty else []
        
    # 3. Configure the metrics to analyze
    metrics_config = [
        ("Dtype Metrics", ['dtype_likely', 'do_strules_dtype', 'should_strules_dtype'], high_dtypes),
        ("Duse Metrics", ['duse_likely', 'do_strules_duse', 'should_strules_duse'], high_duses),
        ("Recipient Metrics", ['recipient_likely', 'do_strules_recipient', 'should_strules_recipient'], high_recipients)
    ]
    
    # --- ASSERTIONS: Ensure all dependent metrics and internal columns exist ---
    for category_name, metrics, high_practices in metrics_config:
        for metric in metrics:
            assert metric in predictor_map, f"Fatal Error: Required metric '{metric}' is missing from predictor_map."
            
            # Check internal dataframe columns
            df_check = predictor_map[metric]
            assert 'topic_condition' in df_check.columns, f"Fatal Error: 'topic_condition' missing in data for '{metric}'."
            assert 'data_practice' in df_check.columns, f"Fatal Error: 'data_practice' missing in data for '{metric}'."

    results = []
    
    # 4. Iterate through combinations and calculate Variance/Stdev
    for app_group_name, app_list in app_groups.items():
        for category_name, metrics, high_practices in metrics_config:
            for metric in metrics:
                df = predictor_map[metric].copy()
                target_col = f"{metric}_target" if f"{metric}_target" in df.columns else metric
                
                df[target_col] = pd.to_numeric(df[target_col], errors='coerce')
                df = df.dropna(subset=[target_col, 'topic_condition', 'data_practice'])
                
                # Filter: Only keep rows belonging to this App Group AND High-Health Practices
                df_filtered = df[(df['topic_condition'].isin(app_list)) & (df['data_practice'].isin(high_practices))]
                
                if not df_filtered.empty:
                    # Method 1: Stdev per practice, then average them
                    stdevs_per_practice = df_filtered.groupby('data_practice')[target_col].std()
                    avg_stdev_across_practices = stdevs_per_practice.mean()
                    
                    # Method 2: Treat all data points as one giant pool and get standard deviation
                    pooled_stdev = df_filtered[target_col].std()
                else:
                    assert False, "No df_filtered found."
                    
                results.append({
                    "App_Healthishness_Level": app_group_name,
                    "Data_Category": category_name,
                    "Metric": metric,
                    "Avg_Stdev_Across_Practices": round(avg_stdev_across_practices, 4) if pd.notna(avg_stdev_across_practices) else "N/A",
                    "Pooled_Stdev": round(pooled_stdev, 4) if pd.notna(pooled_stdev) else "N/A",
                    "N_Qualifying_Practices": df_filtered['data_practice'].nunique()
                })
                
    return pd.DataFrame(results)

# ==========================================
# 8. ORCHESTRATOR HELPER FUNCTIONS
# ==========================================

def extract_hypotheses_from_sheets(sh):
    sheets_map = {
        "v1": sh.worksheet('v1 of Models'),
        "v1_low": sh.worksheet('v1 of Models - Low'),
        "v1_medium": sh.worksheet('v1 of Models - Medium'),
        "v1_high": sh.worksheet('v1 of Models - High'),
        "v1_dummy": sh.worksheet('v1_dummy of Models'),
        "v2": sh.worksheet('v2 of Models'),
        "v3": sh.worksheet('v3 of Models'),
        "v6": sh.worksheet('v6 of Models')
    }

    hypothesis_dict = {}
    for h_id, sheet in sheets_map.items():
        topics = []
        if h_id == "v1_low": topics = ["Flashlight", "Face Editor", "Kinky Dating"]
        elif h_id == "v1_medium": topics = ["STI Dating", "Marijuana Shopping", "Blood Donation"]
        elif h_id == "v1_high": topics = ["Brain Exercises", "Step Counter", "Smart Scale", "STI Help", "Period Tracker"]
        
        print(f"Processing {h_id}...")
        hypothesis_dict[h_id] = extract_hypothesis_special(sheet, h_id, topics)
        
    return hypothesis_dict

def generate_summary_stats(predictor_map, misc_results_dir):
    print("\n" + "="*40)
    print("CALCULATING SUMMARY STATISTICS")
    print("="*40)
    reporter = StatsReporter(predictor_map)
    stats_df = reporter.calculate_stats()
    stats_df.to_csv(misc_results_dir / "summary_statistics_per_data_type.csv", index=False)
    print("Full statistics saved to 'summary_statistics_per_data_type.csv'")

def generate_kruskal_wallis(predictor_map, misc_results_dir):
    print("\n" + "="*40)
    print("CALCULATING PER-APP DIFFERENCES (Kruskal-Wallis)")
    print("="*40)
    tester = StatTester(predictor_map)
    kw_results = tester.check_app_differences(None) 
    if not kw_results.empty:
        kw_results.to_csv(misc_results_dir / "kruskal_wallis_app_differences.csv", index=False)
        print("Kruskal-Wallis results saved to 'kruskal_wallis_app_differences.csv'")

def generate_vif(survey, hypothesis_dict, misc_results_dir):
    print("\n" + "="*40)
    print("CALCULATING VIF FOR EACH HYPOTHESIS")
    print("="*40)
    vif_analyzer = VIFAnalyzer(survey)
    all_vifs = vif_analyzer.analyze_all_hypotheses(hypothesis_dict)
    if not all_vifs.empty:
        all_vifs.to_csv(misc_results_dir / "vif_results_all_hypotheses.csv", index=False)
        print("Full VIF results saved to 'vif_results_all_hypotheses.csv'")

def group_practices_by_prefix(predictor_map):
    """
    Groups the predictor columns and unique evaluated practices into 
    dtype, duse, and recipient lists based on the predictor_map.
    """
    groups = {
        "dtype_columns": set(),
        "duse_columns": set(),
        "recipient_columns": set(),
    }

    for predictor_name, df in predictor_map.items():
        if predictor_name.startswith('dtype_'):
            prefix = "dtype"
        elif predictor_name.startswith('duse_'):
            prefix = "duse"
        elif predictor_name.startswith('recipient_'):
            prefix = "recipient"
        else:
            continue  # Skip any predictors that don't match these prefixes

        # Add the predictor column name to its respective list
        groups[f"{prefix}_columns"].add(predictor_name)
    # Convert sets to sorted lists for clean, predictable output
    return {k: sorted(list(v)) for k, v in groups.items()}

def extract_unique_practices_per_group(predictor_map, groups):
    """
    Extracts all unique 'data_practice' strings for each grouped category 
    (dtype, duse, recipient) from the predictor_map.
    """
    practices_per_group = {}
    
    for group_name, columns in groups.items():
        unique_practices = set()
        
        for col in columns:
            if col in predictor_map:
                # Grab the DataFrame for this specific metric
                df = predictor_map[col]
                
                # Extract all non-null data_practice values and add to our set
                if 'data_practice' in df.columns:
                    practices = df['data_practice'].dropna().unique()
                    unique_practices.update(practices)
            else:
                print(f"Warning: Column '{col}' not found in predictor_map.")
                
        # Convert the set to a sorted list so it's clean and easy to read
        practices_per_group[group_name] = sorted(list(unique_practices))
        
    return practices_per_group

def generate_descriptive_stats(predictor_map, groups, misc_results_dir):
    print("\n" + "="*40)
    print("CALCULATING PER-PRACTICE AVERAGES STRICTLY BY EXTRACTED GROUPS")
    print("="*40)
    
    # 1. Extract the master lists of unique practices for each group
    practice_values_by_group = extract_unique_practices_per_group(predictor_map, groups)
    
    # Map the group keys to their respective health metrics
    group_metric_map = {
        'dtype_columns': {'health': 'dtype_health'},
        'duse_columns': {'health': 'duse_health'},
        'recipient_columns': {'health': 'recipient_health'}
    }

    # --- PER-PRACTICE HEALTHISHNESS STATS ---
    for group_key, metrics in group_metric_map.items():
        health_metric = metrics['health']
        
        if group_key not in practice_values_by_group or health_metric not in predictor_map:
            assert False, f"Fatal Error: {health_metric} - missing from data or groups."
            
        # Get the master list of valid practices for this specific group
        valid_practices = practice_values_by_group[group_key]
        
        # Calculate stats (This already calculates per individual data practice)
        all_stats, high_stats = find_high_average_practices(predictor_map, health_metric, 3.0)
        
        # STRICT FILTER: Only keep rows where the data_practice is in our master list for this group
        all_stats = all_stats[all_stats['data_practice'].isin(valid_practices)].copy()
        high_stats = high_stats[high_stats['data_practice'].isin(valid_practices)].copy()
        
        print(f"\n--- {group_key.upper()} ({health_metric}) ---")
        print(f"Total valid practices evaluated: {len(all_stats)}")
        print(f"Practices with Healthishness >= 3.0: {len(high_stats)}")
        
        assert not high_stats.empty, "No high stats found."
        print(f"\nTop practices for {health_metric}:")
        # Print a preview of the top few practices to the console
        print(high_stats[['data_practice', f'Average_{health_metric}', 'Stdev']].head().to_string(index=False))
        
        # Save the per-practice tables to CSV
        all_stats.to_csv(misc_results_dir / f"per_practice_{health_metric}_all.csv", index=False)
        high_stats.to_csv(misc_results_dir / f"per_practice_{health_metric}_high_only.csv", index=False)

    # --- PER-PRACTICE SENSITIVITY STATS (DTYPE ONLY) ---
    print("\n" + "="*40)
    print("CALCULATING COMPREHENSIVE SENSITIVITY STATS (DTYPE ONLY)")
    print("="*40)
    
    sens_metric = 'dtype_sens'
    group_key = 'dtype_columns'
    
    if group_key in practice_values_by_group and sens_metric in predictor_map:
        valid_practices = practice_values_by_group[group_key]
        
        # Use 0.0 threshold to get everything for dtype_sens
        all_sens_stats, _ = find_high_average_practices(predictor_map, sens_metric, 0.0)
        
        # STRICT FILTER
        all_sens_stats = all_sens_stats[all_sens_stats['data_practice'].isin(valid_practices)].copy()
        
        assert not all_sens_stats.empty, "No all sens stats found."
        print(f"\n--- {group_key.upper()} ({sens_metric}) ---")
        
        # Print a preview showing all the comprehensive stats
        cols_to_print = ['data_practice', f'Average_{sens_metric}', 'Median', 'Mode', 'Stdev', 'Variance']
        existing_cols = [c for c in cols_to_print if c in all_sens_stats.columns]
        print(all_sens_stats[existing_cols].head().to_string(index=False))
        
        all_sens_stats.to_csv(misc_results_dir / f"per_practice_{sens_metric}_all.csv", index=False)
        print(f"\n-> Saved full sensitivity list to 'per_practice_{sens_metric}_all.csv'")
    else:
        assert False, f"Fatal Error: {sens_metric} or {group_key} missing from data."
    
def generate_stdevs(predictor_map, misc_results_dir):
    print("\n" + "="*40)
    print("CALCULATING STANDARD DEVIATIONS")
    print("="*40)
    
    stdev_df = calculate_stdevs(predictor_map)
    
    assert not stdev_df.empty, "No stdevs found."
    print("\nStdev Results (Preview):")
    print(stdev_df.head(10).to_string(index=False))
    stdev_df.to_csv(misc_results_dir / "high_healthishness_stdevs.csv", index=False)
    print("\nFull breakdown saved to 'high_healthishness_stdevs.csv'")

def generate_levenes_test(predictor_map, misc_results_dir):
    levene_results = perform_levenes_test(predictor_map)
    assert not levene_results.empty, "No levene results found."
    print("\nLevene's Test Results (Preview):")
    print(levene_results[['Metric', 'Levene_Statistic', 'p-value', 'Significance']].head(10).to_string(index=False))
    levene_results.to_csv(misc_results_dir / "levenes_test_variance_differences.csv", index=False)
    print("\nLevene's test results saved to 'levenes_test_variance_differences.csv'")

# ==========================================
# 9. MAIN EXECUTION
# ==========================================
def main():
    # 1. Authenticate and Setup Google Sheets
    scopes = ["https://www.googleapis.com/auth/spreadsheets", "https://www.googleapis.com/auth/drive"]
    creds, _ = google.auth.default(scopes=scopes)
    gc = gspread.authorize(creds)

    data_final = gc.open('All Data - Final with Real Practices').worksheet('All Data - Final with Real Practices')
    sh = gc.open('P2S2 Hypotheses')

    # 2. Extract Data
    hypothesis_dict = extract_hypotheses_from_sheets(sh)
    
    print("\nExtracting predictors from raw data...")
    predictor_map = extract_predictors(data_final)

    # 3. Prepare CLMM Data
    print("\nGenerating CLMM Data and Metadata...")
    survey = Survey(predictor_map)
    survey.prepare_clmm_data(hypothesis_dict)
    print("Metadata saved to metadata/metadata.json")

    # 4. Create output directory
    misc_results = Path("misc_results")
    misc_results.mkdir(exist_ok=True)

    groups = group_practices_by_prefix(predictor_map)
    generate_descriptive_stats(predictor_map, groups, misc_results)

    # 5. Run the suite of analyses cleanly
    generate_kruskal_wallis(predictor_map, misc_results)
    generate_vif(survey, hypothesis_dict, misc_results)
    
    generate_stdevs(predictor_map, misc_results)
    generate_levenes_test(predictor_map, misc_results)

    

    print("\nAll Python Pre-Processing Complete.")

if __name__ == '__main__':
    main()