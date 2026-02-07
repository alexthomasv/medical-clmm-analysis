import pandas as pd
import numpy as np

def check_rank_deficiency(file_path):
    # 1. Load the data
    try:
        df = pd.read_csv(file_path)
        print(f"Successfully loaded {len(df)} rows from {file_path}.\n")
    except FileNotFoundError:
        print(f"Error: Could not find file '{file_path}'. Make sure it exists.")
        return

    # List of your fixed effect variables (predictors)
    # These are the ones involved in the "design is column rank deficient" error
    predictors = ['dtype_health', 'purpose_sens', 'dtype_relevant', 
                  'purpose_health', 'dtype_sens']
    
    # Also check against your grouping variable
    group_col = 'topic_condition'

    print("--- DIAGNOSTIC REPORT ---\n")

    # CHECK 1: Constant Columns (Zero Variance)
    # If a column has only 1 unique value, it cannot predict anything.
    print("1. CHECKING FOR CONSTANT COLUMNS:")
    constant_found = False
    for col in predictors:
        if col in df.columns:
            if df[col].nunique() <= 1:
                print(f"   [!] CRITICAL: '{col}' has only 1 unique value: {df[col].unique()}")
                constant_found = True
    
    if not constant_found:
        print("   -> No constant columns found (Good).")

    print("\n" + "-"*30 + "\n")

    # CHECK 2: Perfect Correlation (Multicollinearity)
    # If two columns have correlation > 0.99 (or <-0.99), they are redundant.
    print("2. CHECKING FOR PERFECT CORRELATIONS:")
    
    # Filter only numeric predictors present in the dataframe
    numeric_preds = [c for c in predictors if c in df.columns and pd.api.types.is_numeric_dtype(df[c])]
    
    if len(numeric_preds) > 1:
        corr_matrix = df[numeric_preds].corr().abs()
        
        # Create a mask to ignore the diagonal (self-correlation is always 1)
        mask = np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)
        high_corr = corr_matrix.where(mask).stack()
        
        # Find pairs with correlation > 0.99
        perfect_pairs = high_corr[high_corr > 0.99]
        
        if not perfect_pairs.empty:
            for (col1, col2), val in perfect_pairs.items():
                print(f"   [!] CRITICAL: '{col1}' and '{col2}' are perfectly correlated (r={val:.2f})")
                print(f"       -> You must remove one of these from your model.")
        else:
            print("   -> No perfect correlations found between numeric predictors (Good).")
    else:
        print("   -> Not enough numeric predictors to check correlation.")

    print("\n" + "-"*30 + "\n")

    # CHECK 3: Aliasing with Group (topic_condition)
    # If a predictor is constant within a group, it might conflict with random effects
    # or be redundant if topic_condition was treated as fixed.
    print(f"3. CHECKING ALIASING WITH '{group_col}':")
    
    if group_col in df.columns:
        aliased_found = False
        for col in predictors:
            if col in df.columns:
                # Group by topic and count unique values of the predictor
                counts = df.groupby(group_col)[col].nunique()
                
                # If the maximum unique count is 1, it means for EVERY topic,
                # this variable has only 1 value (it is a property of the topic).
                if counts.max() == 1:
                    print(f"   [!] WARNING: '{col}' is perfectly predicted by '{group_col}'.")
                    print(f"       -> Inside every topic, '{col}' never changes.")
                    aliased_found = True
        
        if not aliased_found:
            print(f"   -> Predictors vary independently of '{group_col}' (Good).")
    else:
        print(f"   -> Column '{group_col}' not found in dataset.")

    print("\n" + "="*30)
    print("SUMMARY ADVICE:")
    print("If you saw [!] CRITICAL errors above, remove those variables from your formula.")
    print("If you saw [!] WARNINGs about aliasing, ensure your random effects are specified correctly.")

# Run the function
check_rank_deficiency('hypothesis_data/v3/v3_0.csv')