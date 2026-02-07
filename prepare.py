# Loads v1/v2 hypothesis

from scipy.stats import kstest, norm

# from google.colab import auth
# auth.authenticate_user()
import gspread
import google.auth

import re
from collections import defaultdict
import pandas as pd


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
    # Remove trailing punctuation and parentheses content like 1207(3635) -> 1207
    # but if you want to keep that info, handle differently.
    # For now, let's just keep the '1207(3635)' as a raw string or parse out 1207.
    # We'll do a minimal approach: remove any parentheses content entirely.
    cleaned = re.sub(r"\(.*?\)", "", string_val)       # remove (3635)
    cleaned = re.sub(r"[^\d.+\-eE]", "", cleaned)      # remove leftover punctuation
    if not cleaned:
        return None

    # Validate if cleaned is still a float
    if not re.match(r"^[+\-]?\d+(\.\d+)?([eE][+\-]?\d+)?$", cleaned):
        return None
    return float(cleaned)

def flatten_summary_dict(nested_dict):
  """
  Flattens a nested dictionary into a single-level dictionary.

  Args:
      nested_dict (dict): A nested dictionary.

  Returns:
      dict: A flattened dictionary with concatenated keys.
  """
  flat_dict = {}

  def flatten(current_key, value):
      if isinstance(value, dict):
          # Recursively flatten the dictionary
          for k, v in value.items():
              new_key = f"{current_key} - {k}" if current_key else k
              flatten(new_key, v)
      else:
          # Assign the value to the flat dictionary
          flat_dict[current_key] = value

  # Start the recursive flattening
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
    # catch variations of "I've never used one"
    elif 'never' in s.lower():
        return 'Never Used One'
    return s  # Fallback (e.g. for empty strings or NO_DATA)

class Hypothesis:
  stmt = ""
  num = 0
  formula = ""
  ind_vars = []
  dep_vars = []
  interactions = []
  def __init__(self, ind_vars, dep_vars, interactions, num, formula):
    self.ind_vars = ind_vars
    self.dep_vars = dep_vars
    self.interactions = interactions
    self.num = num
    self.formula = formula

  def get_base_predictors(self):
      filter = {'usefreq_hf', 'experience_hf', 'gender', 'ethnicity_combined', 'age_bracket', 'region_broad'}
      return set(self.ind_vars) - filter


  def get_r_form(self):
      dep_var_str = self.dep_vars[0]

      # main effects
      ind_parts = list(self.ind_vars)

      # add interactions explicitly
      for inter in self.interactions:
          ind_parts.append('*'.join(inter))

      ind_var_str = ' + '.join(ind_parts) if ind_parts else "1"  # intercept only if empty
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

def extract_hypothesis_special(hypothesis_final, hypothesis_sheet_name):
  hypothesis_list = []
  hypothesis_rows = hypothesis_final.get_all_values()
  h_idx = 0
  for row in hypothesis_rows[1:]:
    # print("row", row)
    # print("row[2]", row[2])
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

    # print("potential_ind_vars", potential_ind_vars)

    for ind_var in potential_ind_vars:
      if '*' in ind_var:
        id = ind_var.split('*')
        assert len(id) == 2 or len(id) == 3, f"{id} {len(id)}"
        id = [i.strip() for i in id]
        interactions.append(tuple(id))
      else:
        ind_vars.append(ind_var)

    # print(ind_vars, dep_vars, formula)
    # print("ind_vars", ind_vars)
    hypothesis_list.append(Hypothesis(ind_vars, dep_vars, interactions, hypothesis_sheet_name + "_" + str(h_idx), formula))
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
    usefreq_hf           = data_row[1494]
    usefreq_med          = data_row[1495]


    experience_hf        = recode_experience(data_row[1496])
    experience_med       = recode_experience(data_row[1497])


    gender               = data_row[1498]
    sex_fwiw             = data_row[1499]
    ethnicity_complex    = data_row[1500]
    ethnicity_simplified = data_row[1501]
    ethnicity_combined   = data_row[1502]
    age                  = data_row[1503]
    age_bracket          = data_row[1504]
    state                = data_row[1505]
    region_granular      = data_row[1506]
    region_broad         = data_row[1507]
    devices_used         = data_row[1508]

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

    # print(data)
    df = pd.DataFrame(data, columns=cols)



    valid_cats = [x for x in ["Goodish", "Baddish", "Never Used One"] if x in df['experience_hf'].unique()]

    df['experience_hf'] = pd.Categorical(
        df['experience_hf'],
        categories=valid_cats,
        ordered=False
    )


    df_sorted = df.sort_values(by=['data_practice',  'topic_condition', 'participant_number'])
    # print(predictor)
    # print(df_sorted)
    predictor_map[predictor] = df_sorted
    # print(df_sorted)
    # assert False

  return predictor_map




class Survey:
  questions = []
  participants = []
  apps = []

  def __init__(self, predictor_map):
    self.predictor_map = predictor_map

  def create_local_table(self, hypothesis: Hypothesis):
    print(hypothesis)
    exceptions = ['purpose_health', 'purpose_sens']
    special = ['region_broad', 'usefreq_med', 'experience_hf', 'experience_med', 'gender', 'ethnicity_combined', 'usefreq_hf', 'age_bracket']


    ind_vars = []
    ind_vars.extend(hypothesis.ind_vars)
    dep_var = hypothesis.dep_vars[0]

    for id in hypothesis.interactions:
      for i in id:
        ind_vars.append(i)

    ind_vars = list(set(ind_vars))

    print(f"Ind vars: {ind_vars}")
    print(f"Dep var: {dep_var}")


    exception_predictors = []
    filtered_predictors = []

    for var in ind_vars:
      if var in special:
        continue
      if var in exceptions:
        exception_predictors.append(var)
      else:
        filtered_predictors.append(var)

    print(filtered_predictors)
    ind_vars_df = [self.predictor_map[ind_var] for ind_var in filtered_predictors]

    # print(ind_vars_df)
    exception_predictors_df = [self.predictor_map[exception].drop('data_practice', axis=1) for exception in exception_predictors]



    dep_var_df = self.predictor_map[dep_var]
    columns = ind_vars_df + [dep_var_df]

    # Reduce to merge all together
    # for c in columns:
    #   print(c)
    final = reduce(lambda left, right: pd.merge(left, right, on=['participant_number', 'topic_condition', 'data_practice'] + special, how='inner'), columns)

    if exception_predictors_df:
      final = reduce(
          lambda L, R: pd.merge(L, R, on=['participant_number','topic_condition'] + special, how='inner'),
          [final] + exception_predictors_df
      )
    print("**** final *****")
    print(final)
    # assert False

    return final

  def run_ks_test(self, hypothesis_list):
    p_value_ref = 0.05
    # First, extract dependent variables
    dependent_variables = [hypothesis.dep_vars[0] for hypothesis in hypothesis_list]
    for var in dependent_variables:
      _, result = self.get(var)
      result = [int(value) for value in result]
      ks_statistic, p_value = kstest(result, 'norm')
      # Perform the one-sample K-S test against a normal distribution
      # Here, norm.cdf specifies the cumulative distribution function of the normal distribution
      print(f"{var}, {ks_statistic}, {p_value}, {p_value < p_value_ref}")


  def prepare_clmm_data(self, hypothesis_dict, do_full):
    from pathlib import Path
    # Get current directory
    cwd = Path.cwd()
    file_prefix_data = cwd / "hypothesis_data"
    file_prefix_metadata = cwd / "metadata"
    import json

    # Create directories if they don't exist
    os.makedirs(file_prefix_data, exist_ok=True)
    os.makedirs(file_prefix_metadata, exist_ok=True)

    metadata = {}
    for hypothesis_name, hypothesis_list in hypothesis_dict.items():
        metadata[hypothesis_name] = {}
        for hypothesis in hypothesis_list:
            hypothesis_name_local = hypothesis.num
            print(f"Hypothesis: {hypothesis_name_local}")
            local_table = self.create_local_table(hypothesis)
            os.makedirs(f"{file_prefix_data}/{hypothesis_name}", exist_ok=True)
            local_table.to_csv(f"{file_prefix_data}/{hypothesis_name}/{hypothesis.num}.csv", index=False)
            # Convert the set returned by get_base_predictors() to a list for JSON serialization
            metadata[hypothesis_name][hypothesis_name_local] = {"ind_vars": hypothesis.ind_vars, "dep_var": hypothesis.dep_vars[0], "interactions": hypothesis.interactions, "base_predictors": list(hypothesis.get_base_predictors()), "do_full": do_full}
            
    # Save to file
    with open(f"{file_prefix_metadata}/metadata.json", "w") as f:
        json.dump(metadata, f, indent=4)


if __name__ == '__main__':
    # This will pick up the credentials from your local environment
    # (set by `gcloud auth application-default login`)
    scopes = [
        "https://www.googleapis.com/auth/spreadsheets",
        "https://www.googleapis.com/auth/drive"
    ]

    # 2. Get credentials with these specific scopes
    creds, _ = google.auth.default(scopes=scopes)
    gc = gspread.authorize(creds)

    # Replace 'your_sheet_name' with the actual name of your Google Sheet
    data_final = gc.open('All Data - Final with Real Practices').worksheet('All Data - Final with Real Practices')

    sh = gc.open('P2S2 Hypotheses')
    hypothesis_final_v1 = sh.worksheet('v1 of Models')
    hypothesis_final_v2 = sh.worksheet('v2 of Models')
    hypothesis_final_v3 = sh.worksheet('v3 of Models')
    hypothesis_final_v1_dummy = sh.worksheet('v1_dummy of Models')
    hypothesis_final_v4 = sh.worksheet('v4 of Models')
    hypothesis_final_v6 = sh.worksheet('v6 of Models')

    hypothesis_id = ["v1", "v1_dummy", "v2", "v3", "v4", "v6"]
    all_hypothesis = [hypothesis_final_v1, hypothesis_final_v1_dummy, hypothesis_final_v2, hypothesis_final_v3, hypothesis_final_v4, hypothesis_final_v6]


    hypothesis_dict = {}
    for hypothesis_name, hypothesis in zip(hypothesis_id, all_hypothesis):
        hypothesis_dict[hypothesis_name] = extract_hypothesis_special(hypothesis, hypothesis_name)

    predictor_map = extract_predictors(data_final)
    survey = Survey(predictor_map)
    # survey.run_ks_test(hypothesis_list)
    survey.prepare_clmm_data(hypothesis_dict, False)
