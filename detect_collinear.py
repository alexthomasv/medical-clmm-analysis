import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from numpy.linalg import cond, svd

def thorough_check(file_path):
    print(f"--- FULL DIAGNOSTIC: {file_path} ---")
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error loading file: {e}")
        return

    # 1. Select the exact predictors involved in the model
    # Update this list if your variable names differ
    predictors = ['dtype_health', 'purpose_sens', 'dtype_relevant', 
                  'purpose_health', 'dtype_sens']
    
    # Check if columns exist
    missing = [c for c in predictors if c not in df.columns]
    if missing:
        print(f"Error: Columns not found: {missing}")
        return

    # Drop NAs (Linear models drop rows with ANY missing data)
    df_clean = df[predictors].dropna()
    X = df_clean.astype(float) # Treat as numeric for the test
    
    print(f"Data Shape: {X.shape[0]} rows x {X.shape[1]} columns\n")

    # --- CHECK 1: CORRELATION MATRIX ---
    print("1. PAIRWISE CORRELATIONS (Threshold > 0.85)")
    corr = X.corr().abs()
    # Mask diagonal
    np.fill_diagonal(corr.values, 0)
    
    high_corr = corr.unstack().sort_values(ascending=False)
    high_corr = high_corr[high_corr > 0.85]
    
    # Remove duplicates (A-B is same as B-A)
    seen = set()
    found_corr = False
    for idx, val in high_corr.items():
        pair = tuple(sorted(idx))
        if pair not in seen:
            seen.add(pair)
            print(f"   [!] WARNING: {pair[0]} <-> {pair[1]} : r = {val:.4f}")
            found_corr = True
    
    if not found_corr:
        print("   -> No dangerously high pairwise correlations found.")

    # --- CHECK 2: VARIANCE INFLATION FACTOR (VIF) ---
    print("\n2. VARIANCE INFLATION FACTOR (VIF)")
    print("   (VIF > 5 is high; VIF > 10 is critical)")
    
    # We add a constant (intercept) because VIF requires it
    X_const = sm.add_constant(X)
    vif_data = pd.DataFrame()
    vif_data["Variable"] = X_const.columns
    vif_data["VIF"] = [variance_inflation_factor(X_const.values, i) 
                       for i in range(X_const.shape[1])]
    
    # Filter out the 'const' row for cleaner output
    print(vif_data[vif_data["Variable"] != "const"].to_string(index=False))

    # --- CHECK 3: EXACT LINEAR COMBINATIONS (R-Squared) ---
    print("\n3. LINEAR COMBINATION CHECK (Predicting X from Others)")
    problem_found = False
    for target in predictors:
        others = [c for c in predictors if c != target]
        y = X[target]
        x_sub = sm.add_constant(X[others])
        
        model = sm.OLS(y, x_sub).fit()
        r2 = model.rsquared
        
        if r2 > 0.99:
            print(f"   [!] CRITICAL: '{target}' is perfectly predicted by others (R2={r2:.4f})")
            problem_found = True
        elif r2 > 0.9:
            print(f"   [!] WARNING: '{target}' is strongly predicted by others (R2={r2:.4f})")
    
    if not problem_found:
        print("   -> No variable is a perfect linear combination of the others.")

    # --- CHECK 4: MATRIX CONDITION NUMBER ---
    print("\n4. GLOBAL STABILITY (Condition Number)")
    # Condition number of the correlation matrix is a standard stability metric
    # < 100 is good. > 1000 is bad.
    c_num = cond(X.corr()) 
    print(f"   Condition Number: {c_num:.2f}")
    
    if c_num < 50:
        print("   -> PASSED. The matrix is extremely stable.")
    elif c_num < 100:
        print("   -> PASSED. The matrix is stable enough.")
    else:
        print("   -> FAILED. The matrix is unstable/singular.")

thorough_check('hypothesis_data/v3/v3_0.csv')