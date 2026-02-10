# Loads v1/v2 hypothesis

from scipy.stats import kstest, norm, kruskal 
# NEW: Import statsmodels for VIF
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant

# from google.colab import auth
# auth.authenticate_user()
import gspread
import google.auth

import re
from collections import defaultdict
import pandas as pd
import numpy as np 

from logging import exception
from functools import reduce
import os

# Define utility code
def get_nickname(string):
  pattern = r"_\d$"  # Matches underscore followed by digit at the end
  extracted_string = re.sub(pattern, "", string)
  return extracted_string

def parse_float_maybe(string_val):
    """
    Tries to parse a float from string_val, stripping non-numeric trailing punctuation.
    Returns None if parsing fails.
    """
    cleaned = re.sub(r"\(.*?\)", "", string_val)       # remove (3635)
    cleaned = re.sub(r"[^\d.+\-eE]", "", cleaned)      # remove leftover punctuation
    if not cleaned:
        return None

    # Validate if cleaned is still a float
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
    """
    Combines 1/2 -> Baddish, 3/4 -> Goodish.
    Standardizes 'Never Used One'.
    """
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

class Hypothesis:
  stmt = ""
  num = 0
  formula = ""
  ind_vars = []
  dep_vars = []
  interactions = []
  def __init__(self, ind_vars, dep_vars, interactions, num, formula, app_topics):
    self.ind_vars = ind_vars
    self.dep_vars = dep_vars
    self.interactions = interactions
    self.num = num
    self.formula = formula
    self.app_topics = app_topics

  def get_base_predictors(self):
      filter = {'usefreq_hf', 'experience_hf', 'gender', 'ethnicity_combined', 'age_bracket', 'region_broad'}
      return set(self.ind_vars) - filter


  def get_r_form(self):
      dep_var_str = self.dep_vars[0]
      ind_parts = list(self.ind_vars)
      for inter in self.interactions:
          ind_parts.append('*'.join(inter))

      ind_var_str = ' + '.join(ind_parts) if ind_parts else "1"
      return f"{dep_var_str} ~ {ind_var_str} ({self.num})"

  def __repr__(self):
    ret_str = self.get_r_form()
    return ret_str

class Predictor:
  predictor = ""
  df = None
  def __init__(self, predictor, df):
    self.predictor = predictor
    self.df = df

def extract_hypothesis_special(hypothesis_final, hypothesis_sheet_name, app_topics):
  hypothesis_list = []
  hypothesis_rows = hypothesis_final.get_all_values()
  h_idx = 0
  for row in hypothesis_rows[1:]:
    ind_vars = row[2].split("+")
    potential_ind_vars = []
    for x in ind_vars:
      x = x.strip()
      if has_paren(x):
        continue
      potential_ind_vars.append(x)

    dep_vars = [row[1].strip()]
    formula = row[2]

    ind_vars = []
    interactions = []

    for ind_var in potential_ind_vars:
      if '*' in ind_var:
        id = ind_var.split('*')
        assert len(id) == 2 or len(id) == 3, f"{id} {len(id)}"
        id = [i.strip() for i in id]
        interactions.append(tuple(id))
      else:
        ind_vars.append(ind_var)

    hypothesis_list.append(Hypothesis(ind_vars, dep_vars, interactions, hypothesis_sheet_name + "_" + str(h_idx), formula, app_topics))
    h_idx += 1

  return hypothesis_list

def extract_predictors(data):
  predictor_map = {}
  data_rows = data.get_all_values()
  predictor_to_data = defaultdict(set)
  predictor_names = data_rows[0][2:]
  data_practices = data_rows[1][2:]
  for data_row in data_rows[2:]:
    participant_number = data_row[0]
    app_topic = data_row[1]
    
    # Safe checks for index out of bounds if rows are short
    if len(data_row) < 1509: 
        continue

    usefreq_hf           = data_row[1494]
    usefreq_med          = data_row[1495]
    experience_hf        = recode_experience(data_row[1496])
    experience_med       = recode_experience(data_row[1497])
    gender               = data_row[1498]
    ethnicity_combined   = recode_ethnicity_combined(data_row[1502])
    age_bracket          = data_row[1504]
    region_broad         = data_row[1507]

    for predictor, data_practice, val in zip(predictor_names, data_practices, data_row[2:]):
      if not val:
        continue
      predictor_to_data[predictor].add((participant_number,
                                        app_topic,
                                        usefreq_hf,
                                        usefreq_med,
                                        experience_hf,
                                        experience_med,
                                        gender,
                                        ethnicity_combined,
                                        age_bracket,
                                        region_broad,
                                        data_practice,
                                        val))


  for predictor in predictor_to_data:
    data = list(predictor_to_data[predictor])
    cols = [
            'participant_number',
            'topic_condition',
            'usefreq_hf',
            'usefreq_med',
            'experience_hf',
            'experience_med',
            'gender',
            'ethnicity_combined',
            'age_bracket',
            'region_broad',
            'data_practice',
            predictor
    ]
    if cols.count(predictor) > 1:
      cols[-1] = f"{predictor}_target"

    df = pd.DataFrame(data, columns=cols)

    valid_cats = [x for x in ["Goodish", "Baddish", "Never Used One"] if x in df['experience_hf'].unique()]

    df['experience_hf'] = pd.Categorical(
        df['experience_hf'],
        categories=valid_cats,
        ordered=False
    )

    df_sorted = df.sort_values(by=['data_practice',  'topic_condition', 'participant_number'])
    predictor_map[predictor] = df_sorted

  return predictor_map

class Survey:
  questions = []
  participants = []
  apps = []

  def __init__(self, predictor_map):
    self.predictor_map = predictor_map

  def create_local_table(self, hypothesis: Hypothesis, topic_conditions_filter):
    # print(hypothesis) # Commented out to reduce noise during bulk VIF run
    exceptions = ['purpose_health', 'purpose_sens']
    special = ['region_broad', 'usefreq_med', 'experience_hf', 'experience_med', 'gender', 'ethnicity_combined', 'usefreq_hf', 'age_bracket']

    ind_vars = []
    ind_vars.extend(hypothesis.ind_vars)
    dep_var = hypothesis.dep_vars[0]

    for id in hypothesis.interactions:
      for i in id:
        ind_vars.append(i)

    ind_vars = list(set(ind_vars))

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
    from pathlib import Path
    cwd = Path.cwd()
    file_prefix_data = cwd / "hypothesis_data"
    file_prefix_metadata = cwd / "metadata"
    import json

    os.makedirs(file_prefix_data, exist_ok=True)
    os.makedirs(file_prefix_metadata, exist_ok=True)

    metadata = {}
    for hypothesis_name, hypothesis_list in hypothesis_dict.items():
        metadata[hypothesis_name] = {}
        for hypothesis in hypothesis_list:
            hypothesis_name_local = hypothesis.num
            # print(f"Hypothesis: {hypothesis_name_local}")
            local_table = self.create_local_table(hypothesis, hypothesis.app_topics)
            os.makedirs(f"{file_prefix_data}/{hypothesis_name}", exist_ok=True)
            local_table.to_csv(f"{file_prefix_data}/{hypothesis_name}/{hypothesis.num}.csv", index=False)
            metadata[hypothesis_name][hypothesis_name_local] = {"ind_vars": hypothesis.ind_vars, "dep_var": hypothesis.dep_vars[0], "interactions": hypothesis.interactions, "base_predictors": list(hypothesis.get_base_predictors())}
            
    with open(f"{file_prefix_metadata}/metadata.json", "w") as f:
        json.dump(metadata, f, indent=4)

# ==========================================
# STATS REPORTER
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
# STATISTICAL TESTER
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
# NEW VIF ANALYZER
# ==========================================
class VIFAnalyzer:
    def __init__(self, survey_instance):
        self.survey = survey_instance

    def compute_vif_for_hypothesis(self, df, ind_vars):
        """
        Calculates VIF for the provided DataFrame and list of independent variables.
        Automatically handles categorical encoding.
        """
        # 1. Filter columns to only those in the hypothesis + available in DF
        # This handles cases where 'ind_vars' might list main effects but not interactions
        valid_predictors = [v for v in ind_vars if v in df.columns]
        
        if not valid_predictors:
            return pd.DataFrame()

        # 2. Subset data
        X = df[valid_predictors].copy()

        # 3. Handle Categoricals: One-Hot Encode (Drop First to avoid dummy trap)
        cat_cols = X.select_dtypes(include=['object', 'category']).columns
        if len(cat_cols) > 0:
            X = pd.get_dummies(X, columns=cat_cols, drop_first=True, dtype=float)

        # 4. Drop missing values
        X = X.dropna()
        if X.empty:
            return pd.DataFrame()

        # 5. Add constant (required for VIF calculation)
        try:
            X = sm.add_constant(X, has_constant='add')
        except Exception:
            # Fallback if add_constant fails for some reason
            X['const'] = 1.0

        # 6. Calculate VIF
        vif_data = []
        for i in range(X.shape[1]):
            feature_name = X.columns[i]
            if feature_name == 'const':
                continue # Skip intercept
            
            try:
                val = variance_inflation_factor(X.values, i)
                vif_data.append({
                    'Predictor': feature_name,
                    'VIF': round(val, 4)
                })
            except Exception as e:
                vif_data.append({
                    'Predictor': feature_name,
                    'VIF': str(e)
                })

        return pd.DataFrame(vif_data)

    def analyze_all_hypotheses(self, hypothesis_dict):
        results = []
        print("\n--- Starting VIF Analysis per Hypothesis ---")
        
        for group_name, h_list in hypothesis_dict.items():
            for h in h_list:
                # Re-create the specific dataset for this hypothesis
                # (Same logic as prepare_clmm_data)
                df = self.survey.create_local_table(h, h.app_topics)
                
                # We calculate VIF based on the independent variables defined in the hypothesis
                # Note: This calculates VIF on main effects. If you have explicit interaction
                # terms (A*B), standard VIF usually looks at the inputs (A, B).
                vif_df = self.compute_vif_for_hypothesis(df, h.ind_vars)
                
                if not vif_df.empty:
                    vif_df['Hypothesis Group'] = group_name
                    vif_df['Hypothesis Num'] = h.num
                    # Reorder cols
                    cols = ['Hypothesis Group', 'Hypothesis Num', 'Predictor', 'VIF']
                    results.append(vif_df[cols])
        
        if results:
            return pd.concat(results, ignore_index=True)
        return pd.DataFrame()


if __name__ == '__main__':
    scopes = [
        "https://www.googleapis.com/auth/spreadsheets",
        "https://www.googleapis.com/auth/drive"
    ]

    creds, _ = google.auth.default(scopes=scopes)
    gc = gspread.authorize(creds)

    data_final = gc.open('All Data - Final with Real Practices').worksheet('All Data - Final with Real Practices')

    sh = gc.open('P2S2 Hypotheses')
    hypothesis_final_v1 = sh.worksheet('v1 of Models')
    
    hypothesis_final_v1_low = sh.worksheet('v1 of Models - Low')
    hypothesis_final_v1_medium = sh.worksheet('v1 of Models - Medium')
    hypothesis_final_v1_high = sh.worksheet('v1 of Models - High')
    
    hypothesis_final_v2 = sh.worksheet('v2 of Models')
    hypothesis_final_v3 = sh.worksheet('v3 of Models')
    hypothesis_final_v1_dummy = sh.worksheet('v1_dummy of Models')
    hypothesis_final_v4 = sh.worksheet('v4 of Models')
    hypothesis_final_v6 = sh.worksheet('v6 of Models')

    hypothesis_id = ["v1", "v1_low", "v1_medium", "v1_high", "v1_dummy", "v2", "v3", "v4", "v6"]
    all_hypothesis = [hypothesis_final_v1, hypothesis_final_v1_low, hypothesis_final_v1_medium, hypothesis_final_v1_high, hypothesis_final_v1_dummy, hypothesis_final_v2, hypothesis_final_v3, hypothesis_final_v4, hypothesis_final_v6]

    hypothesis_dict = {}
    for hypothesis_name, hypothesis in zip(hypothesis_id, all_hypothesis):
        if hypothesis_name == "v1_low":
            app_topics = ["Flashlight", "Face Editor", "Kinky Dating"]
        elif hypothesis_name == "v1_medium":
            app_topics = ["STI Dating", "Marijuana Shopping", "Blood Donation"]
        elif hypothesis_name == "v1_high":
            app_topics = ["Brain Exercises", "Step Counter", "Smart Scale", "STI Help", "Period Tracker"]
        else:
            app_topics = []

        hypothesis_dict[hypothesis_name] = extract_hypothesis_special(hypothesis, hypothesis_name, app_topics)

    predictor_map = extract_predictors(data_final)

    print(predictor_map.keys())
    
    # --- 1. Run Survey Data Prep ---
    survey = Survey(predictor_map)
    survey.prepare_clmm_data(hypothesis_dict)

    # --- 2. Run Statistics Calculation ---
    print("\n" + "="*40)
    print("CALCULATING COMPREHENSIVE STATISTICS")
    print("="*40)
    
    reporter = StatsReporter(predictor_map)
    stats_df = reporter.calculate_stats()
    stats_df.to_csv("summary_statistics_per_data_type.csv", index=False)
    # pd.set_option('display.max_rows', 100) # Optional display
    # print(stats_df.head(30))
    print(f"\nFull statistics saved to 'summary_statistics_per_data_type.csv'")

    # --- 3. Run Kruskal-Wallis Tests ---
    print("\n" + "="*40)
    print("CALCULATING PER-APP DIFFERENCES (Kruskal-Wallis)")
    print("="*40)

    tester = StatTester(predictor_map)
    kw_results = tester.check_app_differences(None) 

    if not kw_results.empty:
        # print(kw_results)
        kw_results.to_csv("kruskal_wallis_app_differences.csv", index=False)
        print("\nKruskal-Wallis results saved to 'kruskal_wallis_app_differences.csv'")
    else:
        print("No results generated.")

    # --- 4. Run VIF Analysis (NEW) ---
    print("\n" + "="*40)
    print("CALCULATING VIF FOR EACH HYPOTHESIS")
    print("="*40)

    vif_analyzer = VIFAnalyzer(survey)
    all_vifs = vif_analyzer.analyze_all_hypotheses(hypothesis_dict)

    if not all_vifs.empty:
        all_vifs.to_csv("vif_results_all_hypotheses.csv", index=False)
        print("\nVIF Analysis Complete.")
        print(all_vifs.head(20)) # Preview
        print("Full VIF results saved to 'vif_results_all_hypotheses.csv'")
    else:
        print("No VIF results generated. Check hypothesis variable names or data availability.")