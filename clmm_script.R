library(jsonlite)
library(ordinal)
library(insight)
library(parallel)

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

format_var <- function(var, data, linear_only = TRUE, linear_method = "numeric") {
  if (!is.ordered(data[[var]])) {
    return(sprintf("`%s`", var))
  }

  if (isTRUE(linear_only)) {
    # --- OPTION A: SIMPLIFY TO LINEAR ---
    if (linear_method == "contrast") {
      return(sprintf("C(`%s`, poly, 1)", var))
    } else {
      # User requested scale(as.numeric(...)) for main effects
      return(sprintf("scale(as.numeric(`%s`))", var))
    }
  } else {
    # --- OPTION B: FULL POLYNOMIAL (L, Q, C) ---
    return(sprintf("`%s`", var))
  }
}

format_as_factor <- function(var, data) {
  sprintf("as.factor(`%s`)", var)
}

build_interaction_terms <- function(interaction_list, data, linear_2way = FALSE, linear_3way = FALSE, linear_method = "numeric") {
  if (is.null(interaction_list) || length(interaction_list) == 0) return(character(0))
  if (!is.list(interaction_list)) return(character(0))
  
  terms <- character(0)
  for (i in seq_along(interaction_list)) {
    components <- as.character(interaction_list[[i]])
    degree <- length(components)
    
    is_linear <- FALSE
    if (degree == 2 && linear_2way) is_linear <- TRUE
    if (degree >= 3 && linear_3way) is_linear <- TRUE
    
    formatted_comps <- sapply(components, function(v) {
      if (!v %in% names(data)) return(v) 
      
      if (!is.ordered(data[[v]])) {
        stop(sprintf("ASSERTION FAILED: Interaction component '%s' is not ordered.", v))
      }

      if (isTRUE(is_linear)) {
         if (linear_method == "contrast") {
            # UPDATED: Use poly, 1 to match format_var
            return(sprintf("C(`%s`, poly, 1)", v))
         } else {
            # User provided code used simple as.numeric for interactions (no scale)
            return(sprintf("scale(as.numeric(`%s`))", v))
         }
      } else {
         return(sprintf("`%s`", v))
      }
    })
    
    terms <- c(terms, paste(formatted_comps, collapse = ":"))
  }
  return(terms)
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

build_data_sub <- function(data, dep_var, all_vars_needed) {
  vars_used <- unique(c(dep_var, all_vars_needed, "topic_condition", "participant_number", "data_practice"))
  vars_used <- vars_used[vars_used %in% names(data)]
  df <- data[, vars_used, drop = FALSE]
  
  bad_never  <- "I've never used one"
  if ("experience_hf"  %in% names(df)) df <- df[df$experience_hf  != bad_never   & !is.na(df$experience_hf), , drop = FALSE]
  if ("experience_med" %in% names(df)) df <- df[df$experience_med != bad_never   & !is.na(df$experience_med), , drop = FALSE]
  
  df[stats::complete.cases(df[, vars_used, drop = FALSE]), , drop = FALSE]
}

flatten_metadata_wrapper <- function(metadata) {
  flat <- list()
  for (group in names(metadata)) {
    if (is.list(metadata[[group]]) && !("dep_var" %in% names(metadata[[group]]))) {
       for (hnum in names(metadata[[group]])) {
         flat[[hnum]] <- metadata[[group]][[hnum]]
       }
    } else {
       flat[[group]] <- metadata[[group]]
    }
  }
  return(flat)
}

get_all_vars_from_hyp <- function(hyp_Obj) {
    base_p <- unlist(hyp_Obj[["base_predictors"]])
    cand_p <- unlist(hyp_Obj[["candidate_terms_predictors"]])
    base_i <- unique(unlist(hyp_Obj[["base_interactions"]]))
    cand_i <- unique(unlist(hyp_Obj[["candidate_terms_interactions"]]))
    unique(c(base_p, cand_p, base_i, cand_i))
}

RANDOM_EFFECTS <- "(1 | topic_condition) + (1 | participant_number) + (1 | data_practice)"

get_data_file_path <- function(hnum, data_prefix) {
  suite_name <- sub("_[0-9]+$", "", hnum)
  file.path(data_prefix, suite_name, paste0(hnum, ".csv"))
}

load_and_prepare_data <- function(data_file_path, dep_var, all_vars) {
  data <- read.csv(data_file_path, check.names = FALSE)
  for (v in all_vars) if (v %in% names(data)) data[[v]] <- cast_by_name(data[[v]], v)
  data[[dep_var]] <- factor(data[[dep_var]], ordered = TRUE)
  data_fixed <- build_data_sub(data, dep_var, all_vars)
  list(raw = data, fixed = data_fixed)
}

apply_poly_contrasts <- function(data_fixed, dep_var) {
  for (v in names(data_fixed)) {
    if (v != dep_var && is.ordered(data_fixed[[v]])) {
      contrasts(data_fixed[[v]]) <- contr.poly(nlevels(data_fixed[[v]]))
    }
  }
  data_fixed
}

make_clmm_formula <- function(dep_var, rhs_str) {
  stats::as.formula(paste(dep_var, "~", rhs_str, "+", RANDOM_EFFECTS))
}

can_add_interaction <- function(candidate_inter, curr_preds, curr_inters) {
    # Check that all main effects are present
    if (!all(candidate_inter %in% curr_preds)) return(FALSE)
    # For 3-way+: check all lower-order sub-interactions are present
    if (length(candidate_inter) > 2) {
        curr_strs <- sapply(curr_inters, function(x) paste(sort(x), collapse=":"))
        for (degree in 2:(length(candidate_inter) - 1)) {
            sub_combs <- combn(candidate_inter, degree, simplify = FALSE)
            for (sub in sub_combs) {
                sub_str <- paste(sort(sub), collapse = ":")
                if (!(sub_str %in% curr_strs)) return(FALSE)
            }
        }
    }
    return(TRUE)
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
    data_file_path <- get_data_file_path(hnum, data_prefix)
    if (!file.exists(data_file_path)) next
    cat(sprintf(">>> Checking Hypothesis: %s\n", hnum), file = result_file, append = TRUE)

    hyp_Obj <- flat_metadata[[hnum]]
    dep_var <- hyp_Obj[["dep_var"]]
    all_vars <- get_all_vars_from_hyp(hyp_Obj)

    prepared <- load_and_prepare_data(data_file_path, dep_var, all_vars)
    data <- prepared$raw; data_fixed <- prepared$fixed
    vars_to_test <- unique(c(unlist(hyp_Obj[["base_predictors"]]), unlist(hyp_Obj[["candidate_terms_predictors"]])))
    
    get_term_str <- function(v) format_var(v, data, linear_only = TRUE, linear_method = "numeric")
    all_terms <- vapply(vars_to_test, get_term_str, character(1))
    rhs_str <- paste(all_terms, collapse = " + ")
    if (rhs_str == "") rhs_str <- "1"
    f_clm_base <- as.formula(paste(dep_var, "~", rhs_str))
    
    fit_clm_base <- tryCatch(ordinal::clm(f_clm_base, data = data_fixed), error = function(e) NULL)
    
    if (is.null(fit_clm_base)) { 
        cat("   [FATAL] Base CLM failed.\n\n", file = result_file, append = TRUE)
        next 
    }
    
    run_single_assumption_test <- function(test_name, test_type) {
        cat(sprintf("   [%s]\n", test_name), file = result_file, append = TRUE)
        for (v in vars_to_test) {
            f_test <- as.formula(paste("~", get_term_str(v)))
            args <- list(f_clm_base, data = data_fixed)
            args[[test_type]] <- f_test
            fit_test <- tryCatch(do.call(ordinal::clm, args), error = function(e) NULL)
            if (is.null(fit_test)) {
                cat(sprintf("       %-20s: [FAIL] Model failed\n", v), file = result_file, append = TRUE)
            } else {
                res <- anova(fit_clm_base, fit_test)
                pval <- res$`Pr(>Chisq)`[2]
                sig <- if (!is.na(pval) && pval < 0.05) "**VIOLATION**" else "Pass"
                cat(sprintf("       %-20s: p=%.4f [%s]\n", v, pval, sig), file = result_file, append = TRUE)
            }
        }
    }

    run_single_assumption_test("A] Proportional Odds Test (Proxy via CLM)", "nominal")
    run_single_assumption_test("B] Scale/Variance Test (Proxy via CLM)", "scale")
    cat("\n", file = result_file, append = TRUE)
  }
}

# ==============================================================================
# 3. PHASE 2: FORWARD SELECTION (HIERARCHY ENFORCED)
# ==============================================================================

run_aic_selection <- function(metadata, data_prefix, target_groups, target_hypotheses, output_root, linear_method = "numeric") {

  json_output_path <- file.path(output_root, "best_models.json")
  if (file.exists(json_output_path)) {
      message(sprintf("INFO: Found existing %s. Skipping Selection.", json_output_path))
      return()
  }

  if (!is.null(target_groups)) { valid_groups <- intersect(names(metadata), target_groups); metadata <- metadata[valid_groups] }
  flat_metadata <- flatten_metadata_wrapper(metadata)

  if (!is.null(target_hypotheses)) {
      valid_hyps <- intersect(names(flat_metadata), target_hypotheses)
      flat_metadata <- flat_metadata[valid_hyps]
      message(sprintf("Filtered to specific hypotheses: %s", paste(valid_hyps, collapse=", ")))
  }

  best_models_map <- list()
  global_log_file <- file.path(output_root, "aic_selection_summary.csv")

  mc <- max(1L, if (is.null(getOption("mc.cores"))) parallel::detectCores(logical=TRUE)-1L else as.integer(getOption("mc.cores")))
  message(sprintf("Running Forward Selection (Linear Method: %s) with %d cores...", linear_method, mc))

  for (hnum in names(flat_metadata)) {
    gc()
    data_file_path <- get_data_file_path(hnum, data_prefix)
    if (!file.exists(data_file_path)) { message(sprintf("Skipping %s (No Data)", hnum)); next }
    message(sprintf("\n>>> Forward selection for: %s", hnum))

    # --- SETUP PER-HYPOTHESIS OUTPUT ---
    hyp_dir <- file.path(output_root, hnum)
    if (!dir.exists(hyp_dir)) dir.create(hyp_dir, recursive = TRUE)
    history_file <- file.path(hyp_dir, "aic_history.csv")
    verbose_file <- file.path(hyp_dir, "aic_history_verbose.txt")

    cat("Step,Action,AIC,Status,Candidates_Tested,Best_Candidate_AIC\n", file = history_file)
    cat(sprintf("VERBOSE FORWARD SELECTION LOG FOR %s\n", hnum), file = verbose_file)
    cat("===================================================\n", file = verbose_file, append=TRUE)

    hyp_Obj <- flat_metadata[[hnum]]
    dep_var <- hyp_Obj[["dep_var"]]
    base_preds <- unlist(hyp_Obj[["base_predictors"]])
    base_inters <- hyp_Obj[["base_interactions"]]
    cand_preds_pool <- unlist(hyp_Obj[["candidate_terms_predictors"]])
    cand_inters_pool <- hyp_Obj[["candidate_terms_interactions"]]

    all_vars <- get_all_vars_from_hyp(hyp_Obj)
    prepared <- load_and_prepare_data(data_file_path, dep_var, all_vars)
    data <- prepared$raw; data_fixed <- prepared$fixed

    if (linear_method == "contrast") data_fixed <- apply_poly_contrasts(data_fixed, dep_var)

    build_rhs <- function(kept_preds, kept_inters_list) {
        main_terms <- vapply(kept_preds, function(v) format_var(v, data, linear_only=TRUE, linear_method=linear_method), character(1))
        inter_terms <- build_interaction_terms(kept_inters_list, data, linear_2way=TRUE, linear_3way=TRUE, linear_method=linear_method)
        rhs_vec <- c(main_terms, inter_terms)
        if (length(rhs_vec) == 0) return("1")
        paste(rhs_vec, collapse = " + ")
    }

    control <- ordinal::clmm.control(gradTol = 1e-6)

    # -----------------------------------------------------------
    # 1. FIT BASE MODEL
    # -----------------------------------------------------------
    curr_preds <- base_preds
    curr_inters <- base_inters
    rhs_str <- build_rhs(curr_preds, curr_inters)
    f_base <- make_clmm_formula(dep_var, rhs_str)

    fit <- tryCatch(ordinal::clmm(f_base, data = data_fixed, Hess = FALSE, method = "nlminb", control = control), error = function(e) NULL)

    base_ok <- FALSE
    if (!is.null(fit)) {
        if (!any(abs(coef(fit)) > 10)) {
            base_ok <- TRUE
        }
    }

    if (!base_ok) {
        message(sprintf("  -> Base model failed or unstable for %s. Skipping.", hnum))
        cat(sprintf("0,Base Model,NA,FAILED,0,NA\n"), file = history_file, append = TRUE)
        cat(sprintf("Base model FAILED for %s. Skipping hypothesis.\n", hnum), file = verbose_file, append = TRUE)
        next
    }

    best_aic <- AIC(fit)
    cat(sprintf("0,Base Model,%.2f,Accepted,0,NA\n", best_aic), file = history_file, append = TRUE)
    cat(sprintf("\nStep 0: Base Model | AIC: %.2f\n", best_aic), file = verbose_file, append = TRUE)
    message(sprintf("  -> Base AIC: %.2f", best_aic))

    # -----------------------------------------------------------
    # 2. FORWARD SELECTION LOOP
    # -----------------------------------------------------------
    step_count <- 1
    improved <- TRUE

    while (improved) {
      improved <- FALSE

      # Collect eligible candidates: remaining predictors + hierarchy-valid interactions
      avail_preds <- setdiff(cand_preds_pool, curr_preds)

      curr_inter_strs <- sapply(curr_inters, function(x) paste(sort(x), collapse=":"))
      valid_int_indices <- integer(0)
      for (i in seq_along(cand_inters_pool)) {
          candidate <- cand_inters_pool[[i]]
          cand_str <- paste(sort(candidate), collapse=":")
          if (cand_str %in% curr_inter_strs) next
          if (can_add_interaction(candidate, curr_preds, curr_inters)) {
              valid_int_indices <- c(valid_int_indices, i)
          }
      }

      # Build job list
      jobs <- list()
      for (v in avail_preds) jobs[[length(jobs) + 1]] <- list(type = "var", val = v)
      for (idx in valid_int_indices) jobs[[length(jobs) + 1]] <- list(type = "intr", val = cand_inters_pool[[idx]])

      if (length(jobs) == 0) {
          cat(sprintf("%d,Stop (No Candidates),%.2f,Final,0,NA\n", step_count, best_aic), file = history_file, append = TRUE)
          cat(sprintf("\nStep %d: No eligible candidates remain. Stopping.\n", step_count), file = verbose_file, append = TRUE)
          break
      }

      cat(sprintf("\nStep %d: Testing %d candidates (current AIC: %.2f)\n", step_count, length(jobs), best_aic), file = verbose_file, append = TRUE)

      # Fit all candidates in parallel
      trial_results <- parallel::mclapply(jobs, function(job) {
          try_preds <- curr_preds; try_inters <- curr_inters
          if (job$type == "var") { try_preds <- c(try_preds, job$val) } else { try_inters <- c(try_inters, list(job$val)) }
          rhs_try <- build_rhs(try_preds, try_inters)
          f_try <- make_clmm_formula(dep_var, rhs_try)

          fit_try <- tryCatch(ordinal::clmm(f_try, data = data_fixed, Hess = FALSE, method = "nlminb", control = control), error = function(e) NULL)

          if (is.null(fit_try)) return(list(job = job, status = "ERROR", aic = Inf))

          # Safety check: reject models with exploding coefficients (separation)
          if (any(abs(coef(fit_try)) > 10)) {
              return(list(job = job, status = "REJECTED (Separation/Exploding Coef)", aic = Inf))
          }

          list(job = job, status = "OK", aic = AIC(fit_try), preds = try_preds, inters = try_inters)
      }, mc.cores = mc)

      # Log all candidate results to verbose file
      for (res in trial_results) {
          if (is.null(res)) next
          cand_name <- if (res$job$type == "var") res$job$val else paste(res$job$val, collapse = ":")
          aic_str <- if (is.infinite(res$aic)) "NA" else sprintf("%.2f", res$aic)
          cat(sprintf("  Candidate: %-30s | Status: %-40s | AIC: %s\n", cand_name, res$status, aic_str), file = verbose_file, append = TRUE)
      }

      # Filter to valid results and find best
      valid_results <- Filter(function(x) !is.null(x) && x$status == "OK" && is.finite(x$aic), trial_results)

      candidates_tested_str <- paste(sapply(jobs, function(j) if (j$type == "var") j$val else paste(j$val, collapse = ":")), collapse = "; ")

      if (length(valid_results) == 0) {
          cat(sprintf("%d,Stop (All Candidates Failed),%.2f,Final,%d,NA\n", step_count, best_aic, length(jobs)), file = history_file, append = TRUE)
          cat(sprintf("  All %d candidates failed or were rejected.\n", length(jobs)), file = verbose_file, append = TRUE)
          break
      }

      aics <- vapply(valid_results, function(x) x$aic, numeric(1))
      best_idx <- which.min(aics)
      trial_best_aic <- aics[best_idx]
      best_try <- valid_results[[best_idx]]
      added_name <- if (best_try$job$type == "var") best_try$job$val else paste(best_try$job$val, collapse = ":")

      if (trial_best_aic < best_aic) {
          improvement <- best_aic - trial_best_aic
          curr_preds <- best_try$preds
          curr_inters <- best_try$inters
          best_aic <- trial_best_aic
          message(sprintf("  -> Step %d: Added %s (AIC: %.2f)", step_count, added_name, best_aic))
          cat(sprintf("%d,Added %s,%.2f,Accepted,%d,%.2f\n", step_count, added_name, best_aic, length(jobs), trial_best_aic), file = history_file, append = TRUE)
          cat(sprintf("  >>> ACCEPTED: %s (AIC: %.2f, improvement: %.2f)\n", added_name, trial_best_aic, improvement), file = verbose_file, append = TRUE)
          improved <- TRUE
      } else {
          cat(sprintf("%d,Stop (No Improvement),%.2f,Final,%d,%.2f\n", step_count, best_aic, length(jobs), trial_best_aic), file = history_file, append = TRUE)
          cat(sprintf("  No improvement. Best candidate: %s (AIC: %.2f vs current %.2f). Stopping.\n", added_name, trial_best_aic, best_aic), file = verbose_file, append = TRUE)
      }
      step_count <- step_count + 1
    }

    # -----------------------------------------------------------
    # 3. RECORD FINAL MODEL
    # -----------------------------------------------------------
    message(sprintf("  -> Final model AIC: %.2f | Predictors: %s", best_aic, paste(curr_preds, collapse = ", ")))
    cat(sprintf("\n===================================================\nFINAL MODEL: AIC=%.2f\nPredictors: %s\nInteractions: %s\n",
        best_aic, paste(curr_preds, collapse = ", "),
        if (length(curr_inters) > 0) paste(sapply(curr_inters, function(x) paste(x, collapse = ":")), collapse = ", ") else "None"),
        file = verbose_file, append = TRUE)

    best_models_map[[hnum]] <- list(kept_vars = curr_preds, kept_intrs_list = curr_inters)
    intr_strs <- sapply(curr_inters, function(x) paste(x, collapse=":"))
    row_data <- data.frame(hypothesis = hnum, AIC = as.numeric(best_aic), kept_predictors = paste(curr_preds, collapse = "+"), interactions = paste(intr_strs, collapse = "+"))
    write.table(row_data, file = global_log_file, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(global_log_file))
  }

  write_json(best_models_map, json_output_path, pretty = TRUE, auto_unbox = TRUE)
}

# ==============================================================================
# 4. PHASE 3: FINAL CLMM FITTING
# ==============================================================================

perform_linearity_check <- function(final_fit, data, kept_vars, dep_var, control, output_dir) {
  log_path <- file.path(output_dir, "linearity_check.txt")
  cat(sprintf("Linearity Check Report for %s\n", dep_var), file = log_path)
  cat("========================================\n", file = log_path, append = TRUE)
  
  for (v in kept_vars) {
    if (is.ordered(data[[v]])) {
      f_check <- make_clmm_formula(dep_var, format_as_factor(v, data))
      fit_check <- tryCatch(ordinal::clmm(f_check, data = data, Hess = FALSE, method = "nlminb", control = control), error = function(e) NULL)
      if (!is.null(fit_check)) {
        coefs <- coef(fit_check)
        idx <- grep(paste0("as.factor\\(`?", v, "`?\\)"), names(coefs))
        if (length(idx) > 1) {
          r <- cor(1:(length(idx)+1), c(0, coefs[idx]))
          status <- if(r>0.9) "[PASS]" else if(r>0.7) "[WARN]" else "[FAIL]"
          msg <- sprintf("Linearity Check %-15s: %s (r=%.3f)", v, status, r)
          message(msg)
          cat(msg, "\n", file = log_path, append = TRUE)
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
        inter_strs <- trimws(strsplit(df$interactions[i], "\\+")[[1]]); inter_strs <- inter_strs[nzchar(inter_strs)]
        kept_intrs_list <- lapply(inter_strs, function(s) { parts <- strsplit(s, ":")[[1]]; gsub("`", "", parts) })
        best_models_map[[hnum]] <- list(kept_vars = kept_vars, kept_intrs_list = kept_intrs_list)
    }
    return(best_models_map)
}

run_final_fitting <- function(metadata, data_prefix, json_dir, run_dir, 
                              force_linear = FALSE, linear_2way = FALSE, linear_3way = FALSE, linear_method = "numeric") {
    
    # SAFETY CHECK: PRINT CONFIGURATION
    message(sprintf(">>> STARTING FINAL FIT WITH CONFIG: linear_method = '%s'", linear_method))
    
    json_path <- file.path(json_dir, "best_models.json")
    csv_path  <- file.path(json_dir, "aic_selection_summary.csv")
    selected_models <- NULL
    
    if (file.exists(json_path)) {
        message(sprintf("Loading model configs from JSON: %s", json_path))
        selected_models <- read_json(json_path)
    } else if (file.exists(csv_path)) {
        message(sprintf("JSON missing. Reconstructing from CSV: %s", csv_path))
        selected_models <- reconstruct_best_models_from_csv(csv_path)
    }
    
    if (is.null(selected_models)) stop("Error: Could not find 'best_models.json' or 'aic_selection_summary.csv'.")
    meta_flat <- flatten_metadata_wrapper(metadata) 
    valid_models <- intersect(names(selected_models), names(meta_flat))
    message(sprintf("\nRunning Final Fitting for %d models (filtered)...", length(valid_models)))
    
    main_effects_linear <- force_linear
    use_linear_2way     <- if(force_linear) TRUE else linear_2way
    use_linear_3way     <- if(force_linear) TRUE else linear_3way
    
    if(force_linear) {
        message(sprintf(">>> MODE: ALL LINEAR (Method: %s)", linear_method))
    } else {
        message(sprintf(">>> MODE: HYBRID\n    Main Effects: FACTORS\n    2-Way Linear: %s\n    3-Way+ Linear: %s\n    Method: %s",
                use_linear_2way, use_linear_3way, linear_method))
    }
    
    global_stats_file <- file.path(run_dir, "all_model_stats.csv")
    global_coefs_file <- file.path(run_dir, "all_model_coefficients.csv")
    
    for (hnum in valid_models) {
        gc()
        message(sprintf("Fitting: %s", hnum))
        
        hyp_dir <- file.path(run_dir, hnum)
        if (!dir.exists(hyp_dir)) dir.create(hyp_dir, recursive = TRUE)
        log_file <- file.path(hyp_dir, "clmm_output.txt")
        
        model_cfg <- selected_models[[hnum]]
        kept_vars <- as.character(model_cfg$kept_vars)
        kept_intrs_list <- model_cfg$kept_intrs_list 
        
        hyp_Obj <- meta_flat[[hnum]]
        dep_var <- hyp_Obj[["dep_var"]]
        
        data_file_path <- get_data_file_path(hnum, data_prefix)
        all_vars <- unique(c(dep_var, kept_vars, unlist(kept_intrs_list)))
        prepared <- load_and_prepare_data(data_file_path, dep_var, all_vars)
        data <- prepared$raw; data_fixed <- prepared$fixed

        if (linear_method == "contrast") data_fixed <- apply_poly_contrasts(data_fixed, dep_var)
        
        main_terms <- vapply(kept_vars, function(v) format_var(v, data, linear_only = main_effects_linear, linear_method = linear_method), character(1))
        inter_terms <- build_interaction_terms(kept_intrs_list, data, linear_2way = use_linear_2way, linear_3way = use_linear_3way, linear_method = linear_method)
        rhs_vec <- c(main_terms, inter_terms)
        if (length(rhs_vec) == 0) rhs_vec <- "1"
        rhs_str <- paste(rhs_vec, collapse = " + ")
        
        f_final <- make_clmm_formula(dep_var, rhs_str)
        
        # Try with gradTol=1e-6 first, fall back to default if Hessian fails
        control_tight <- ordinal::clmm.control(gradTol = 1e-6)
        fit <- tryCatch(ordinal::clmm(f_final, data = data_fixed, Hess = TRUE, method = "nlminb", control = control_tight), error = function(e) NULL)

        if (!is.null(fit) && any(is.nan(summary(fit)$coefficients[, "Std. Error"]))) {
            message("  -> Hessian failed with gradTol=1e-6, retrying with default control...")
            fit <- tryCatch(ordinal::clmm(f_final, data = data_fixed, Hess = TRUE, method = "nlminb"), error = function(e) NULL)
        }

        if (is.null(fit)) { message("  -> Final fit failed."); next }
        message(sprintf("  -> Converged. AIC: %.2f", AIC(fit)))
        
        sink(log_file)
        print(summary(fit))
        sink()
        
        row_data <- data.frame(hypothesis = hnum, AIC = as.numeric(AIC(fit)), nobs = fit$info$nobs, max_grad = fit$info$max.grad, formula = paste(deparse(f_final), collapse = ""))
        write.table(row_data, file = global_stats_file, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(global_stats_file))
        
        coef_summary <- summary(fit)$coefficients
        estimates <- coef_summary[, "Estimate"]
        std_errors <- coef_summary[, "Std. Error"]
        det_df <- data.frame(Hypothesis = hnum, Term = rownames(coef_summary), Type = ifelse(grepl("\\|", rownames(coef_summary)), "Threshold", "Predictor"), Estimate = estimates, Std_Error = std_errors, Z_Value = coef_summary[, "z value"], P_Value = coef_summary[, "Pr(>|z|)"], Odds_Ratio = exp(estimates), CI_Lower_95 = exp(estimates - 1.96 * std_errors), CI_Upper_95 = exp(estimates + 1.96 * std_errors), stringsAsFactors = FALSE)
        write.table(det_df, file = global_coefs_file, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(global_coefs_file))
        
        perform_linearity_check(fit, data_fixed, kept_vars, dep_var, control, hyp_dir)
    }
}

# ==============================================================================
# 5. DRIVER
# ==============================================================================

if (sys.nframe() == 0) {  # Only run when executed directly, not when sourced

# --- CONFIGURATION ---
DO_ASSUMPTIONS <- TRUE
DO_SELECTION   <- FALSE
DO_FINAL_FIT   <- TRUE

# Modeling Strategy
FORCE_LINEAR_FINAL_FIT <- TRUE
LINEAR_2WAY_INTERACTIONS <- TRUE   # Force 2-way interactions to be linear
LINEAR_3WAY_INTERACTIONS <- TRUE   # Force 3-way interactions to be linear
LINEAR_METHOD <- "contrast"         # Set to "numeric" for 3rd option (as.numeric) or "contrast" for C(x, poly, 1)

# Filtering
TARGET_GROUPS  <- NULL  # Filter by Suite (e.g., "v6")
TARGET_HYPOTHESES <- NULL # Filter by Specific Hypothesis (e.g., "v6_8"). Set NULL for all.
TARGET_PREDICTORS <- NULL 

# --- SETUP ---
file_prefix <- "metadata/"
metadata    <- fromJSON(paste0(file_prefix, "metadata.json"))
data_prefix <- "hypothesis_data/"
timestamp   <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
RESULTS_ROOT <- "results"
run_dir     <- file.path(RESULTS_ROOT, paste0("run_", timestamp))
if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE)

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
      if (any(unlist(hyp$base_predictors) %in% TARGET_PREDICTORS)) kept_hypotheses[[h]] <- hyp
    }
    if (length(kept_hypotheses) > 0) filtered_metadata[[g]] <- kept_hypotheses
  }
  metadata <- filtered_metadata
  message(sprintf("  -> Kept models in %d groups.", length(metadata)))
  if(length(metadata) == 0) stop("No models found matching criteria.")
}

# --- SOURCE CHECK LOGIC: find latest best_models.json across all runs ---
json_source_dir <- NULL
existing_jsons <- Sys.glob(file.path(RESULTS_ROOT, "run_*", "best_models.json"))
if (length(existing_jsons) > 0) {
    # run_ dirs are timestamped so lexicographic sort gives chronological order
    json_source_dir <- dirname(sort(existing_jsons, decreasing = TRUE)[1])
    message(sprintf("Found existing 'best_models.json' in %s.", json_source_dir))
    if (DO_SELECTION) {
        message("Skipping Step 2 (Model Selection) and using existing configuration.")
        DO_SELECTION <- FALSE
    }
} else {
    message("No existing 'best_models.json' found in any run. Will run selection from scratch (if enabled).")
}

# --- PIPELINE ---
if (DO_ASSUMPTIONS) {
  message("--- Step 1: Assumption Tests ---")
  run_all_assumption_tests(metadata, data_prefix, NULL, run_dir)
}

if (DO_SELECTION) {
  message("\n--- Step 2: Model Selection (Forward AIC / Linear) ---")
  # UPDATED: Pass LINEAR_METHOD config
  run_aic_selection(metadata, data_prefix, TARGET_GROUPS, TARGET_HYPOTHESES, run_dir, linear_method = LINEAR_METHOD)
}

if (DO_FINAL_FIT) {
  message("\n--- Step 3: Final Fitting (CLMM) ---")
  # Use run_dir if selection just ran, otherwise the latest run with best_models.json
  source_dir <- if (DO_SELECTION) run_dir else json_source_dir
  if (is.null(source_dir)) stop("Cannot run final fitting: no best_models.json found in any run directory.")
  run_final_fitting(metadata, data_prefix, source_dir, run_dir,
                    force_linear = FORCE_LINEAR_FINAL_FIT,
                    linear_2way = LINEAR_2WAY_INTERACTIONS,
                    linear_3way = LINEAR_3WAY_INTERACTIONS,
                    linear_method = LINEAR_METHOD)
}

message("\nAll tasks complete. Results in: ", run_dir)

}  # end if (sys.nframe() == 0)