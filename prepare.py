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
    if len(data_row) < 1509: continue

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

    for predictor, data_practice, val in zip(predictor_names, data_practices, data_row[2:]):
      if not val: continue
      predictor_to_data[predictor].add((participant_number, app_topic, usefreq_hf, usefreq_med, experience_hf, experience_med, gender, ethnicity_combined, age_bracket, region_broad, data_practice, val))

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

  def create_local_table(self, hypothesis: Hypothesis, topic_conditions_filter):
    exceptions = ['purpose_health', 'purpose_sens']
    special = ['region_broad', 'usefreq_med', 'experience_hf', 'experience_med', 'gender', 'ethnicity_combined', 'usefreq_hf', 'age_bracket']
    
    # Identify all needed columns from ind_vars (which includes base + candidates)
    # Also ensure constituents of interactions are present
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
    columns = ind_vars_df + [dep_var_df]

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

        X = df[valid_predictors].copy()
        cat_cols = X.select_dtypes(include=['object', 'category']).columns
        if len(cat_cols) > 0:
            X = pd.get_dummies(X, columns=cat_cols, drop_first=True, dtype=float)

        X = X.dropna()
        if X.empty: return pd.DataFrame()

        try:
            X = sm.add_constant(X, has_constant='add')
        except Exception:
            X['const'] = 1.0

        vif_data = []
        for i in range(X.shape[1]):
            feature_name = X.columns[i]
            if feature_name == 'const': continue
            try:
                val = variance_inflation_factor(X.values, i)
                vif_data.append({'Predictor': feature_name, 'VIF': round(val, 4)})
            except Exception as e:
                vif_data.append({'Predictor': feature_name, 'VIF': str(e)})

        return pd.DataFrame(vif_data)

    def analyze_all_hypotheses(self, hypothesis_dict):
        results = []
        print("\n--- Starting VIF Analysis per Hypothesis ---")
        
        for group_name, h_list in hypothesis_dict.items():
            for h in h_list:
                df = self.survey.create_local_table(h, h.app_topics)
                # Check VIF for Base Predictors (Candidates usually excluded from initial VIF check)
                # Combining base + candidates for VIF check might be too aggressive if they aren't all in the model
                # Let's check BASE predictors by default
                vars_to_check = h.base_predictors 
                if not vars_to_check: vars_to_check = h.ind_vars # Fallback if base is empty
                
                vif_df = self.compute_vif_for_hypothesis(df, vars_to_check)
                
                if not vif_df.empty:
                    vif_df['Hypothesis Group'] = group_name
                    vif_df['Hypothesis Num'] = h.num
                    cols = ['Hypothesis Group', 'Hypothesis Num', 'Predictor', 'VIF']
                    results.append(vif_df[cols])
        
        if results: return pd.concat(results, ignore_index=True)
        return pd.DataFrame()

# ==========================================
# 7. MAIN EXECUTION
# ==========================================
if __name__ == '__main__':
    scopes = ["https://www.googleapis.com/auth/spreadsheets", "https://www.googleapis.com/auth/drive"]
    creds, _ = google.auth.default(scopes=scopes)
    gc = gspread.authorize(creds)

    data_final = gc.open('All Data - Final with Real Practices').worksheet('All Data - Final with Real Practices')
    sh = gc.open('P2S2 Hypotheses')

    # Define Sheets
    sheets_map = {
        "v1": sh.worksheet('v1 of Models'),
        "v1_low": sh.worksheet('v1 of Models - Low'),
        "v1_medium": sh.worksheet('v1 of Models - Medium'),
        "v1_high": sh.worksheet('v1 of Models - High'),
        "v1_dummy": sh.worksheet('v1_dummy of Models'),
        "v2": sh.worksheet('v2 of Models'),
        "v3": sh.worksheet('v3 of Models'),
        # "v4": sh.worksheet('v4 of Models'),
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

    print("\nExtracting predictors from raw data...")
    predictor_map = extract_predictors(data_final)

    print("\nGenerating CLMM Data and Metadata...")
    survey = Survey(predictor_map)
    survey.prepare_clmm_data(hypothesis_dict)
    print("Metadata saved to metadata/metadata.json")

    # --- Run Statistics ---
    print("\n" + "="*40)
    print("CALCULATING SUMMARY STATISTICS")
    print("="*40)
    reporter = StatsReporter(predictor_map)
    stats_df = reporter.calculate_stats()
    stats_df.to_csv("summary_statistics_per_data_type.csv", index=False)
    print("Full statistics saved to 'summary_statistics_per_data_type.csv'")

    # --- Run Kruskal-Wallis ---
    print("\n" + "="*40)
    print("CALCULATING PER-APP DIFFERENCES (Kruskal-Wallis)")
    print("="*40)
    tester = StatTester(predictor_map)
    kw_results = tester.check_app_differences(None) 
    if not kw_results.empty:
        kw_results.to_csv("kruskal_wallis_app_differences.csv", index=False)
        print("Kruskal-Wallis results saved to 'kruskal_wallis_app_differences.csv'")

    # --- Run VIF ---
    print("\n" + "="*40)
    print("CALCULATING VIF FOR EACH HYPOTHESIS")
    print("="*40)
    vif_analyzer = VIFAnalyzer(survey)
    all_vifs = vif_analyzer.analyze_all_hypotheses(hypothesis_dict)
    if not all_vifs.empty:
        all_vifs.to_csv("vif_results_all_hypotheses.csv", index=False)
        print("Full VIF results saved to 'vif_results_all_hypotheses.csv'")
    
    print("\nAll Python Pre-Processing Complete.")