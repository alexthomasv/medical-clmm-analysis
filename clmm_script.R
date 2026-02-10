# install.packages("ordinal")
# install.packages("insight")
library(jsonlite)
library(ordinal)
library(insight)
library(parallel)

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

format_var <- function(var, data, linear_only = TRUE, linear_method = "numeric") {
  # If variable is ordered factor and we want linear constraint:
  if (is.ordered(data[[var]]) && isTRUE(linear_only)) {
    if (linear_method == "contrast") {
      # Method B: Orthogonal Polynomial Linear Contrast
      return(sprintf("C(`%s`, poly, 1)", var))
    } else {
      # Method A: Standard Numeric Scaling (Default)
      return(sprintf("scale(as.numeric(`%s`))", var))
    }
  } 
  # Otherwise return variable name (Factor default)
  else {
    return(sprintf("`%s`", var))
  }
}

format_as_factor <- function(var, data) {
  sprintf("as.factor(`%s`)", var)
}

# UPDATED: Now accepts linear_method to control formatting
build_interactions_dynamic <- function(interactions, data, linear_2way = FALSE, linear_3way = FALSE, linear_method = "numeric") {
  if (is.null(interactions)) return(character(0))
  
  if (is.list(interactions)) {
    if (length(interactions) == 0L) return(character(0))
    lst <- lapply(interactions, function(x) as.character(unlist(x, use.names = FALSE)))
    max_len <- max(lengths(lst))
    mat <- matrix(NA_character_, nrow = length(lst), ncol = max_len)
    for (i in seq_along(lst)) {
      vec <- lst[[i]]
      mat[i, seq_along(vec)] <- vec
    }
    interactions <- mat
  }
  if (is.data.frame(interactions)) interactions <- as.matrix(interactions)
  if (!is.matrix(interactions) || any(dim(interactions) == 0L)) return(character(0))

  vars <- unique(as.vector(interactions))
  vars <- vars[!is.na(vars)]
  vars <- trimws(as.character(vars))
  vars <- vars[nzchar(vars)]
  missing <- setdiff(vars, names(data))
  if (length(missing)) stop("Missing variables in `data`: ", paste(missing, collapse = ", "))

  for (v in vars) {
    if (!is.factor(data[[v]])) data[[v]] <- factor(data[[v]], ordered = TRUE)
  }

  out_list <- apply(interactions, 1, function(r) {
    r <- trimws(as.character(r))
    r <- r[!is.na(r) & nzchar(r)]
    if (!length(r)) return(NULL)
    if (any(duplicated(r))) return(NULL)
    
    degree <- length(r)
    is_linear <- FALSE
    if (degree == 2 && linear_2way) is_linear <- TRUE
    if (degree >= 3 && linear_3way) is_linear <- TRUE
    
    term_pieces <- sapply(r, function(v) {
      if (is.ordered(data[[v]]) && isTRUE(is_linear)) {
         if (linear_method == "contrast") {
            sprintf("C(`%s`, poly, 1)", v)
         } else {
            sprintf("scale(as.numeric(`%s`))", v)
         }
      } else {
        sprintf("`%s`", v)
      }
    })
    list(term_str = paste(term_pieces, collapse = ":"), parents  = r)
  })
  out_list[!sapply(out_list, is.null)]
}

# Legacy wrapper for AIC selection (defaults to all linear, numeric method)
build_pairwise_interactions <- function(interactions, data, linear_only = TRUE) {
    build_interactions_dynamic(interactions, data, linear_2way = linear_only, linear_3way = linear_only, linear_method = "numeric")
}

cast_by_name <- function(x, name) {
  lev_usefreq     <- c("Less than once a year","A few times a year", "A few times a month","A few times a week","Daily")
  lev_gender      <- c("Woman","Man","Non-binary")
  lev_ethnicity   <- c("White","Asian","Black","Latino/Hispanic","Mixed","Other")
  lev_age         <- c("18-20","21-44","45-64","65+")
  lev_region      <- c("South","Northeast","West","Midwest")
  lev_experience_hf <- c("Goodish", "Baddish", "Never Used One")
  lev_experience_med <- c("Goodish", "Baddish", "Never Used One")

  x <- trimws(as.character(x))
  if(name %in% c("gender", "ethnicity_combined", "age_bracket", "region_broad")) {
     x[x == "NO_DATA"] <- NA
  }

  switch(name,
         "usefreq_hf"        = factor(x, levels = lev_usefreq,    ordered = TRUE),
         "usefreq_med"       = factor(x, levels = lev_usefreq,    ordered = TRUE),
         "experience_hf"     = factor(x, levels = lev_experience_hf, ordered = FALSE),
         "experience_med"    = factor(x, levels = lev_experience_med, ordered = FALSE),
         "age_bracket"       = factor(x, levels = lev_age,        ordered = TRUE),
         "gender"            = factor(x, levels = lev_gender,     ordered = FALSE),
         "ethnicity_combined"= factor(x, levels = lev_ethnicity,  ordered = FALSE),
         "region_broad"      = factor(x, levels = lev_region,     ordered = FALSE),
         factor(x, ordered = TRUE)
  )
}

build_data_sub <- function(data, dep_var, kept_vars) {
  vars_used <- unique(c(dep_var, kept_vars, "topic_condition", "participant_number", "data_practice"))
  df <- data[, vars_used, drop = FALSE]
  bad_never  <- "I've never used one"
  if ("experience_hf"  %in% names(df)) df <- df[df$experience_hf  != bad_never   & !is.na(df$experience_hf), , drop = FALSE]
  if ("experience_med" %in% names(df)) df <- df[df$experience_med != bad_never   & !is.na(df$experience_med), , drop = FALSE]
  df[stats::complete.cases(df[, vars_used, drop = FALSE]), , drop = FALSE]
}

# --- HELPER: Flatten Metadata ---
flatten_metadata_wrapper <- function(metadata) {
  flat <- list()
  for (group in names(metadata)) {
    if (is.list(metadata[[group]]) && !("ind_vars" %in% names(metadata[[group]]))) {
       for (hnum in names(metadata[[group]])) {
         flat[[hnum]] <- metadata[[group]][[hnum]]
       }
    } else {
       flat[[group]] <- metadata[[group]]
    }
  }
  return(flat)
}

# ==============================================================================
# 2. PHASE 1: ASSUMPTION TESTING
# ==============================================================================

run_all_assumption_tests <- function(metadata, data_prefix, target_groups, output_dir) {
  result_file <- file.path(output_dir, "assumption_results_detailed.txt")
  cat("=======================================================\n", file = result_file)
  cat(" ASSUMPTION CHECK REPORT (Proxy Tests via CLM)\n", file = result_file, append = TRUE)
  cat("=======================================================\n\n", file = result_file, append = TRUE)

  if (!is.null(target_groups)) { 
      valid_groups <- intersect(names(metadata), target_groups)
      metadata <- metadata[valid_groups]
  }
  flat_metadata <- flatten_metadata_wrapper(metadata)
  
  for (hnum in names(flat_metadata)) {
    suite_name <- sub("_[0-9]+$", "", hnum)
    data_file_path <- file.path(data_prefix, suite_name, paste0(hnum, ".csv"))
    if (!file.exists(data_file_path)) next
    cat(sprintf(">>> Checking Hypothesis: %s\n", hnum), file = result_file, append = TRUE)
    
    hyp_Obj <- flat_metadata[[hnum]]
    ind_vars <- hyp_Obj[["ind_vars"]]; dep_var <- hyp_Obj[["dep_var"]]; base_preds <- hyp_Obj[["base_predictors"]]
    all_vars <- unique(c(ind_vars, base_preds))
    
    data <- read.csv(data_file_path, check.names = FALSE)
    for (v in all_vars) if (v %in% names(data)) data[[v]] <- cast_by_name(data[[v]], v)
    data[[dep_var]] <- factor(data[[dep_var]], ordered = TRUE)
    
    do_linear_only <- !grepl("_do_full$", hnum)
    data_fixed <- build_data_sub(data, dep_var, all_vars)
    vars_to_test <- all_vars
    # Note: Assumption tests always default to numeric scaling for stability in proxy checks
    get_term_str <- function(v) format_var(v, data, linear_only = do_linear_only, linear_method = "numeric")
    
    all_terms <- vapply(vars_to_test, get_term_str, character(1))
    rhs_str <- paste(all_terms, collapse = " + ")
    f_clm_base <- as.formula(paste(dep_var, "~", rhs_str))
    
    fit_clm_base <- tryCatch(ordinal::clm(f_clm_base, data = data_fixed), error = function(e) NULL)
    
    if (is.null(fit_clm_base)) { cat("   [FATAL] Base CLM failed.\n\n", file = result_file, append = TRUE); next }
    cat("   [A] Proportional Odds Test (Proxy via CLM)\n", file = result_file, append = TRUE)
    for (v in vars_to_test) {
        f_nom <- as.formula(paste("~", get_term_str(v)))
        fit_nom <- tryCatch(ordinal::clm(f_clm_base, nominal = f_nom, data = data_fixed), error = function(e) NULL)
        if (is.null(fit_nom)) { cat(sprintf("       %-20s: [FAIL] Model failed\n", v), file = result_file, append = TRUE) } else {
             res <- anova(fit_clm_base, fit_nom); pval <- res$`Pr(>Chisq)`[2]; sig <- if (!is.na(pval) && pval < 0.05) "**VIOLATION**" else "Pass"
             cat(sprintf("       %-20s: p=%.4f [%s]\n", v, pval, sig), file = result_file, append = TRUE)
        }
    }
    cat("   [B] Scale/Variance Test (Proxy via CLM)\n", file = result_file, append = TRUE)
    for (v in vars_to_test) {
        f_scale <- as.formula(paste("~", get_term_str(v)))
        fit_scale <- tryCatch(ordinal::clm(f_clm_base, scale = f_scale, data = data_fixed), error = function(e) NULL)
        if (is.null(fit_scale)) { cat(sprintf("       %-20s: [FAIL] Model failed\n", v), file = result_file, append = TRUE) } else {
             res <- anova(fit_clm_base, fit_scale); pval <- res$`Pr(>Chisq)`[2]; sig <- if (!is.na(pval) && pval < 0.05) "**VIOLATION**" else "Pass"
             cat(sprintf("       %-20s: p=%.4f [%s]\n", v, pval, sig), file = result_file, append = TRUE)
        }
    }
    cat("\n", file = result_file, append = TRUE)
  }
}

# ==============================================================================
# 3. PHASE 2: AIC MODEL SELECTION
# ==============================================================================

run_aic_selection <- function(metadata, data_prefix, target_groups, output_dir) {
  
  json_output_path <- file.path(output_dir, "best_models.json")
  if (file.exists(json_output_path)) {
      message(sprintf("INFO: Found existing %s. Skipping AIC Selection.", json_output_path))
      return()
  }

  if (!is.null(target_groups)) { valid_groups <- intersect(names(metadata), target_groups); metadata <- metadata[valid_groups] }
  flat_metadata <- flatten_metadata_wrapper(metadata)
  
  best_models_map <- list()
  selection_log_file <- file.path(output_dir, "aic_selection_log.csv")
  mc <- max(1L, if (is.null(getOption("mc.cores"))) parallel::detectCores(logical=TRUE)-1L else as.integer(getOption("mc.cores")))
  message(sprintf("Running AIC Selection (Linear Assumption) with %d cores...", mc))

  for (hnum in names(flat_metadata)) {
    gc()
    suite_name <- sub("_[0-9]+$", "", hnum)
    data_file_path <- file.path(data_prefix, suite_name, paste0(hnum, ".csv"))
    if (!file.exists(data_file_path)) { message(sprintf("Skipping %s (No Data)", hnum)); next }
    message(sprintf("Selecting: %s", hnum))
    
    hyp_Obj <- flat_metadata[[hnum]]
    ind_vars <- hyp_Obj[["ind_vars"]]; dep_var <- hyp_Obj[["dep_var"]]; inters <- hyp_Obj[["interactions"]]; base_predictors <- hyp_Obj[["base_predictors"]]
    if (is.null(base_predictors)) base_predictors <- character(0)
    
    all_vars <- unique(c(ind_vars, base_predictors))
    data <- read.csv(data_file_path, check.names = FALSE)
    for (v in all_vars) if (v %in% names(data)) data[[v]] <- cast_by_name(data[[v]], v)
    data[[dep_var]] <- factor(data[[dep_var]], ordered = TRUE)

    # Note: AIC Selection defaults to numeric/linear because it's faster and cleaner for selection
    do_linear_only <- TRUE 
    interaction_defs <- build_pairwise_interactions(inters, data, linear_only = do_linear_only)
    data_fixed <- build_data_sub(data, dep_var, all_vars)

    build_formula <- function(kept_vars, kept_intrs_strs) {
      main_terms <- vapply(kept_vars, function(v) format_var(v, data, linear_only = do_linear_only, linear_method = "numeric"), character(1))
      rhs_terms <- c(main_terms, kept_intrs_strs); rhs_str <- paste(rhs_terms, collapse = " + ")
      stats::as.formula(paste(dep_var, "~", if (nzchar(rhs_str)) rhs_str else "1", "+ (1 | topic_condition) + (1 | participant_number) + (1 | data_practice)"))
    }

    kept_vars  <- base_predictors; kept_intrs <- character(0)
    control <- ordinal::clmm.control(gradTol = 1e-6)
    fit <- tryCatch(ordinal::clmm(build_formula(kept_vars, kept_intrs), data = data_fixed, Hess = FALSE, method = "nlminb", control = control), error = function(e) NULL)
    if (is.null(fit)) { message(sprintf("  -> Base model failed for %s", hnum)); next }
    best_aic <- AIC(fit)

    improved <- TRUE
    while (improved) {
      improved <- FALSE
      cand_vars <- setdiff(all_vars, kept_vars)
      cand_intrs <- list()
      if (length(interaction_defs) > 0) { for (item in interaction_defs) if (!item$term_str %in% kept_intrs && all(item$parents %in% kept_vars)) cand_intrs <- c(cand_intrs, list(item)) }
      jobs <- list()
      for (v in cand_vars) jobs[[length(jobs)+1]] <- list(type="var", val=v)
      for (i in cand_intrs) jobs[[length(jobs)+1]] <- list(type="intr", val=i$term_str)
      if (length(jobs) == 0) break
      
      trial_results <- parallel::mclapply(jobs, function(job) {
          try_vars <- kept_vars; try_intrs <- kept_intrs
          if (job$type == "var") try_vars <- c(try_vars, job$val) else try_intrs <- c(try_intrs, job$val)
          fit_try <- tryCatch(ordinal::clmm(build_formula(try_vars, try_intrs), data = data_fixed, Hess = FALSE, method = "nlminb", control = control), error = function(e) NULL)
          if (is.null(fit_try)) return(NULL)
          list(job=job, aic=AIC(fit_try), vars=try_vars, intrs=try_intrs)
        }, mc.cores = mc)
      trial_results <- Filter(Negate(is.null), trial_results)
      if (length(trial_results) == 0) break
      aics <- vapply(trial_results, function(x) x$aic, numeric(1))
      best_idx <- which.min(aics)
      if (aics[best_idx] < best_aic) {
        best_try <- trial_results[[best_idx]]; kept_vars <- best_try$vars; kept_intrs <- best_try$intrs; best_aic <- best_try$aic
        fit <- ordinal::clmm(build_formula(kept_vars, kept_intrs), data = data_fixed, Hess = FALSE, method = "nlminb", control = control); improved <- TRUE
      }
    }
    best_models_map[[hnum]] <- list(kept_vars = kept_vars, kept_intrs = kept_intrs, interaction_parents = lapply(kept_intrs, function(istr) { for(def in interaction_defs) { if(def$term_str == istr) return(def$parents) }; return(NULL) }))
    row_data <- data.frame(hypothesis = hnum, AIC = as.numeric(best_aic), kept_predictors = paste(kept_vars, collapse = " + "), interactions = paste(kept_intrs, collapse = " + "))
    write.table(row_data, file = selection_log_file, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(selection_log_file))
  }
  write_json(best_models_map, json_output_path, pretty = TRUE, auto_unbox = TRUE)
}

# ==============================================================================
# 4. PHASE 3: FINAL CLMM FITTING
# ==============================================================================

perform_linearity_check <- function(final_fit, data, kept_vars, kept_intrs, dep_var, control) {
  if (length(kept_vars) == 0) return()
  for (v in kept_vars) {
    if (is.ordered(data[[v]])) {
      f_check <- stats::as.formula(paste(dep_var, "~", format_as_factor(v, data), "+ (1 | topic_condition) + (1 | participant_number) + (1 | data_practice)"))
      fit_check <- tryCatch(ordinal::clmm(f_check, data = data, Hess = FALSE, method = "nlminb", control = control), error = function(e) NULL)
      if (!is.null(fit_check)) {
        coefs <- coef(fit_check)
        idx <- grep(paste0("as.factor\\(`?", v, "`?\\)"), names(coefs))
        if (length(idx) > 1) {
          r <- cor(1:(length(idx)+1), c(0, coefs[idx])); status <- if(r>0.9) "[PASS]" else if(r>0.7) "[WARN]" else "[FAIL]"
          message(sprintf("Linearity Check %-15s: %s (r=%.3f)", v, status, r))
        }
      }
    }
  }
}

reconstruct_best_models_from_csv <- function(csv_path) {
    if (!file.exists(csv_path)) return(NULL)
    df <- read.csv(csv_path, stringsAsFactors=FALSE)
    df <- df[!duplicated(df$hypothesis, fromLast=TRUE), ]
    best_models_map <- list()
    for (i in 1:nrow(df)) {
        hnum <- df$hypothesis[i]
        kept_vars <- trimws(strsplit(df$kept_predictors[i], "\\+")[[1]]); kept_vars <- kept_vars[nzchar(kept_vars)]
        kept_intrs <- trimws(strsplit(df$interactions[i], "\\+")[[1]]); kept_intrs <- kept_intrs[nzchar(kept_intrs)]
        parents_list <- lapply(kept_intrs, function(s) {
            parts <- strsplit(s, ":")[[1]]
            clean_parts <- gsub("scale\\(as\\.numeric\\(`?|`?\\)\\)", "", parts)
            clean_parts <- gsub("`", "", clean_parts)
            clean_parts <- gsub("C\\(`?", "", clean_parts) # Remove C( wrapper from new method
            clean_parts <- gsub("`, poly, 1\\)", "", clean_parts)
            return(clean_parts)
        })
        best_models_map[[hnum]] <- list(kept_vars = kept_vars, kept_intrs = kept_intrs, interaction_parents = parents_list)
    }
    return(best_models_map)
}

run_final_fitting <- function(metadata, data_prefix, json_dir, output_dir, 
                              force_linear = FALSE, 
                              linear_2way = FALSE, 
                              linear_3way = FALSE,
                              linear_method = "numeric") {
    
    # 1. Load Selected Models
    json_path <- file.path(json_dir, "best_models.json")
    csv_path  <- file.path(json_dir, "aic_selection_log.csv")
    selected_models <- NULL
    
    if (file.exists(json_path)) {
        message(sprintf("Loading model configs from JSON: %s", json_path))
        selected_models <- read_json(json_path)
    } else if (file.exists(csv_path)) {
        message(sprintf("JSON missing. Reconstructing from CSV: %s", csv_path))
        selected_models <- reconstruct_best_models_from_csv(csv_path)
    }
    
    if (is.null(selected_models)) stop("Error: Could not find 'best_models.json' or 'aic_selection_log.csv'.")
    
    # 2. Setup Data Filter
    meta_flat <- flatten_metadata_wrapper(metadata) 
    valid_models <- intersect(names(selected_models), names(meta_flat))
    
    message(sprintf("\nRunning Final Fitting for %d models (filtered)...", length(valid_models)))
    
    # Configuration Logic
    main_effects_linear <- force_linear
    use_linear_2way     <- if(force_linear) TRUE else linear_2way
    use_linear_3way     <- if(force_linear) TRUE else linear_3way
    
    if(force_linear) {
        message(sprintf(">>> MODE: ALL LINEAR (Main + All Interactions)\n    Method: %s", linear_method))
    } else {
        message(sprintf(">>> MODE: HYBRID\n    Main Effects: %s\n    2-Way Linear: %s\n    3-Way+ Linear: %s\n    Method: %s",
                "FACTORS", use_linear_2way, use_linear_3way, linear_method))
    }
    
    for (hnum in valid_models) {
        gc()
        message(sprintf("Fitting: %s", hnum))
        
        # --- PATHS ---
        hyp_dir <- file.path(output_dir, hnum)
        if (!dir.exists(hyp_dir)) dir.create(hyp_dir, recursive = TRUE)
        summary_file <- file.path(hyp_dir, "model_stats.csv")
        detailed_file <- file.path(hyp_dir, "coefficients.csv")
        log_file <- file.path(hyp_dir, "clmm_output.txt")
        
        # --- LOAD & PREP ---
        model_cfg <- selected_models[[hnum]]
        kept_vars <- as.character(model_cfg$kept_vars)
        kept_intrs_parents <- model_cfg$interaction_parents 
        hyp_Obj <- meta_flat[[hnum]]
        dep_var <- hyp_Obj[["dep_var"]]
        
        suite_name <- sub("_[0-9]+$", "", hnum)
        data_file_path <- file.path(data_prefix, suite_name, paste0(hnum, ".csv"))
        data <- read.csv(data_file_path, check.names = FALSE)
        
        all_vars <- unique(c(dep_var, kept_vars, unlist(kept_intrs_parents)))
        for (v in all_vars) if (v %in% names(data)) data[[v]] <- cast_by_name(data[[v]], v)
        data[[dep_var]] <- factor(data[[dep_var]], ordered = TRUE)
        data_fixed <- build_data_sub(data, dep_var, all_vars)
        
        # --- BUILD FORMULA ---
        # 1. Main Effects
        main_terms <- vapply(kept_vars, function(v) format_var(v, data, linear_only = main_effects_linear, linear_method = linear_method), character(1))
        
        # 2. Interactions (Dynamic based on degree)
        final_intrs_strs <- character(0)
        if (length(kept_intrs_parents) > 0) {
             dummy_inters <- list(); idx<-1; for(p in kept_intrs_parents) { dummy_inters[[idx]]<-p; idx<-idx+1 }
             
             # Call new dynamic builder
             rebuilt_defs <- build_interactions_dynamic(dummy_inters, data, linear_2way = use_linear_2way, linear_3way = use_linear_3way, linear_method = linear_method)
             for(def in rebuilt_defs) final_intrs_strs <- c(final_intrs_strs, def$term_str)
        }
        
        rhs_terms <- c(main_terms, final_intrs_strs)
        rhs_str   <- paste(rhs_terms, collapse = " + ")
        f_final   <- stats::as.formula(paste(dep_var, "~", if (nzchar(rhs_str)) rhs_str else "1", "+ (1 | topic_condition) + (1 | participant_number) + (1 | data_practice)"))
        
        # --- FIT ---
        control <- ordinal::clmm.control(gradTol = 1e-6)
        fit <- tryCatch(ordinal::clmm(f_final, data = data_fixed, Hess = TRUE, method = "nlminb", control = control), error = function(e) NULL)
        
        if (is.null(fit)) { message("  -> Final fit failed."); next }
        message(sprintf("  -> Converged. AIC: %.2f", AIC(fit)))
        
        # --- SAVE OUTPUTS ---
        sink(log_file)
        print(summary(fit))
        sink()
        
        row_data <- data.frame(hypothesis = hnum, AIC = as.numeric(AIC(fit)), nobs = fit$info$nobs, max_grad = fit$info$max.grad, formula = paste(deparse(f_final), collapse = ""))
        write.table(row_data, file = summary_file, append = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
        
        coef_summary <- summary(fit)$coefficients
        estimates <- coef_summary[, "Estimate"]
        std_errors <- coef_summary[, "Std. Error"]
        
        det_df <- data.frame(
          Hypothesis  = hnum,
          Term        = rownames(coef_summary),
          Type        = ifelse(grepl("\\|", rownames(coef_summary)), "Threshold", "Predictor"),
          Estimate    = estimates,
          Std_Error   = std_errors,
          Z_Value     = coef_summary[, "z value"],
          P_Value     = coef_summary[, "Pr(>|z|)"],
          Odds_Ratio  = exp(estimates),
          CI_Lower_95 = exp(estimates - 1.96 * std_errors),
          CI_Upper_95 = exp(estimates + 1.96 * std_errors),
          stringsAsFactors = FALSE
        )
        write.table(det_df, file = detailed_file, append = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)
        
        perform_linearity_check(fit, data_fixed, kept_vars, final_intrs_strs, dep_var, control)
    }
}

# ==============================================================================
# 5. DRIVER
# ==============================================================================

# --- CONFIGURATION ---
DO_ASSUMPTIONS <- FALSE
DO_SELECTION   <- FALSE
DO_FINAL_FIT   <- TRUE

# 1. MODELING STRATEGY
# Force all predictors (Main + Interaction) to be numeric/linear?
FORCE_LINEAR_FINAL_FIT <- FALSE

# Hybrid Toggles (Only active if FORCE_LINEAR_FINAL_FIT is FALSE)
LINEAR_2WAY_INTERACTIONS <- TRUE   # Force 2-way interactions to be linear
LINEAR_3WAY_INTERACTIONS <- TRUE   # Force 3-way interactions to be linear

# Linear Method: "numeric" = scale(as.numeric(x)) | "contrast" = C(x, poly, 1)
LINEAR_METHOD <- "contrast"

# 2. FILTERING
# Group Filter (e.g. c("v3") or NULL to run all)
TARGET_GROUPS  <- c("v6")

# Predictor Filter (e.g. c("purpose_health") or NULL to run all)
TARGET_PREDICTORS <- NULL 

# --- SETUP ---
file_prefix <- "metadata/"
metadata    <- fromJSON(paste0(file_prefix, "metadata.json"))
data_prefix <- "hypothesis_data/"
timestamp   <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
run_dir     <- file.path("results", paste0("run_", timestamp))
if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE)

RESULTS_ROOT <- "results"

# --- FILTERING LOGIC ---
if (!is.null(TARGET_GROUPS)) {
  valid_groups <- intersect(names(metadata), TARGET_GROUPS)
  metadata <- metadata[valid_groups]
  message(sprintf("Filtered to %d groups.", length(metadata)))
}

if (!is.null(TARGET_PREDICTORS)) {
  message(sprintf("Filtering for predictors: %s", paste(TARGET_PREDICTORS, collapse=", ")))
  filtered_metadata <- list()
  for (g in names(metadata)) {
    group_content <- metadata[[g]]
    kept_hypotheses <- list()
    for (h in names(group_content)) {
      hyp <- group_content[[h]]
      if (any(hyp$ind_vars %in% TARGET_PREDICTORS)) kept_hypotheses[[h]] <- hyp
    }
    if (length(kept_hypotheses) > 0) filtered_metadata[[g]] <- kept_hypotheses
  }
  metadata <- filtered_metadata
  message(sprintf("  -> Kept models in %d groups.", length(metadata)))
  if(length(metadata) == 0) stop("No models found matching criteria.")
}

# --- SOURCE CHECK LOGIC ---
json_source_dir <- NULL
if (file.exists(file.path(RESULTS_ROOT, "best_models.json"))) {
    message(sprintf("Found existing 'best_models.json' in %s.", RESULTS_ROOT))
    if(DO_SELECTION) {
        message("Skipping Step 2 (Model Selection) and using existing configuration.")
        DO_SELECTION <- FALSE
    }
    json_source_dir <- RESULTS_ROOT
} else {
    message("No existing 'best_models.json' found. Will run selection from scratch (if enabled).")
    json_source_dir <- run_dir 
}

# --- PIPELINE ---
if (DO_ASSUMPTIONS) {
  message("--- Step 1: Assumption Tests ---")
  run_all_assumption_tests(metadata, data_prefix, NULL, run_dir)
}

if (DO_SELECTION) {
  message("\n--- Step 2: Model Selection (Forward AIC / Linear) ---")
  run_aic_selection(metadata, data_prefix, NULL, run_dir)
  json_source_dir <- run_dir
}

if (DO_FINAL_FIT) {
  message("\n--- Step 3: Final Fitting (CLMM) ---")
  run_final_fitting(metadata, data_prefix, json_source_dir, run_dir, 
                    force_linear = FORCE_LINEAR_FINAL_FIT, 
                    linear_2way = LINEAR_2WAY_INTERACTIONS,
                    linear_3way = LINEAR_3WAY_INTERACTIONS,
                    linear_method = LINEAR_METHOD)
}

message("\nAll tasks complete. Results in: ", run_dir)