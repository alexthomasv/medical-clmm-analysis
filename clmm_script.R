# install.packages("ordinal")
# install.packages("insight")
library(jsonlite)
library(ordinal)
library(insight)
library(parallel)

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

format_var <- function(var, data, linear_only = TRUE) {
  if (is.ordered(data[[var]]) && isTRUE(linear_only)) {
    sprintf("scale(as.numeric(`%s`))", var) 
  } else {
    sprintf("`%s`", var) 
  }
}

format_as_factor <- function(var, data) {
  sprintf("as.factor(`%s`)", var)
}

build_pairwise_interactions <- function(interactions, data, linear_only = TRUE) {
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
    term_pieces <- sapply(r, function(v) {
      if (is.ordered(data[[v]]) && isTRUE(linear_only)) {
        sprintf("scale(as.numeric(`%s`))", v)
      } else {
        sprintf("`%s`", v)
      }
    })
    list(term_str = paste(term_pieces, collapse = ":"), parents  = r)
  })
  out_list[!sapply(out_list, is.null)]
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

# ==============================================================================
# 2. PHASE 1: ASSUMPTION TESTING (Proxy Tests using CLM)
# ==============================================================================

run_all_assumption_tests <- function(metadata, data_prefix, target_groups, output_dir) {
  
  result_file <- file.path(output_dir, "assumption_results_detailed.txt")
  cat("=======================================================\n", file = result_file)
  cat(" ASSUMPTION CHECK REPORT (Proxy Tests via CLM)\n", file = result_file, append = TRUE)
  cat(" Note: Random effects are dropped for these tests.\n", file = result_file, append = TRUE)
  cat("=======================================================\n\n", file = result_file, append = TRUE)

  # Flatten metadata
  if (!is.null(target_groups)) {
    valid_groups <- intersect(names(metadata), target_groups)
    metadata <- metadata[valid_groups]
  }
  
  flat_metadata <- list()
  for (group in names(metadata)) {
    if (is.list(metadata[[group]]) && !("ind_vars" %in% names(metadata[[group]]))) {
       for (hnum in names(metadata[[group]])) flat_metadata[[hnum]] <- metadata[[group]][[hnum]]
    } else {
       flat_metadata[[group]] <- metadata[[group]]
    }
  }
  metadata <- flat_metadata
  
  for (hnum in names(metadata)) {
    suite_name <- sub("_[0-9]+$", "", hnum)
    data_file_path <- file.path(data_prefix, suite_name, paste0(hnum, ".csv"))
    
    if (!file.exists(data_file_path)) next
    
    cat(sprintf(">>> Checking Hypothesis: %s\n", hnum), file = result_file, append = TRUE)
    
    # 1. Load & Prep Data
    ind_vars <- metadata[[hnum]][["ind_vars"]]
    dep_var  <- metadata[[hnum]][["dep_var"]]
    base_preds <- metadata[[hnum]][["base_predictors"]]
    all_vars <- unique(c(ind_vars, base_preds))
    
    data <- read.csv(data_file_path, check.names = FALSE)
    for (v in all_vars) if (v %in% names(data)) data[[v]] <- cast_by_name(data[[v]], v)
    data[[dep_var]] <- factor(data[[dep_var]], ordered = TRUE)
    
    do_linear_only <- !grepl("_do_full$", hnum)
    data_fixed <- build_data_sub(data, dep_var, all_vars)
    vars_to_test <- all_vars
    
    get_term_str <- function(v) format_var(v, data, linear_only = do_linear_only)
    all_terms <- vapply(vars_to_test, get_term_str, character(1))
    rhs_str <- paste(all_terms, collapse = " + ")
    
    # ---------------------------------------------------------
    # SETUP: Base CLM Model (No Random Effects)
    # ---------------------------------------------------------
    # We use CLM because CLMM does not support 'nominal' or 'scale' arguments
    f_clm_base <- as.formula(paste(dep_var, "~", rhs_str))
    
    fit_clm_base <- tryCatch(ordinal::clm(f_clm_base, data = data_fixed), error = function(e) paste("Error:", conditionMessage(e)))
    
    if (is.character(fit_clm_base)) {
      cat(sprintf("   [FATAL] Base CLM failed: %s\n\n", fit_clm_base), file = result_file, append = TRUE)
      next
    }

    # ---------------------------------------------------------
    # TEST A: Nominal Test (Proportional Odds)
    # ---------------------------------------------------------
    cat("   [A] Proportional Odds Test (Proxy via CLM)\n", file = result_file, append = TRUE)
    
    for (v in vars_to_test) {
        # Formula for Nominal: ~ predictor
        f_nom <- as.formula(paste("~", get_term_str(v)))
        
        fit_nom <- tryCatch(
          ordinal::clm(f_clm_base, nominal = f_nom, data = data_fixed), 
          error = function(e) { return(list(status="error", msg=conditionMessage(e))) },
          warning = function(w) { return(list(status="warning", msg=conditionMessage(w))) }
        )
        
        if (is.list(fit_nom) && !is.null(fit_nom$status)) {
             cat(sprintf("       %-20s: [FAIL] %s\n", v, fit_nom$msg), file = result_file, append = TRUE)
        } else {
             res <- anova(fit_clm_base, fit_nom)
             pval <- res$`Pr(>Chisq)`[2]
             sig <- if (!is.na(pval) && pval < 0.05) "**VIOLATION**" else "Pass"
             cat(sprintf("       %-20s: p=%.4f [%s]\n", v, pval, sig), file = result_file, append = TRUE)
        }
    }

    # ---------------------------------------------------------
    # TEST B: Scale Test (Heteroscedasticity)
    # ---------------------------------------------------------
    cat("   [B] Scale/Variance Test (Proxy via CLM)\n", file = result_file, append = TRUE)
    
    for (v in vars_to_test) {
        # Formula for Scale: ~ predictor
        f_scale <- as.formula(paste("~", get_term_str(v)))
        
        fit_scale <- tryCatch(
          ordinal::clm(f_clm_base, scale = f_scale, data = data_fixed), 
          error = function(e) { return(list(status="error", msg=conditionMessage(e))) },
          warning = function(w) { return(list(status="warning", msg=conditionMessage(w))) }
        )
        
        if (is.list(fit_scale) && !is.null(fit_scale$status)) {
             cat(sprintf("       %-20s: [FAIL] %s\n", v, fit_scale$msg), file = result_file, append = TRUE)
        } else {
             res <- anova(fit_clm_base, fit_scale)
             pval <- res$`Pr(>Chisq)`[2]
             sig <- if (!is.na(pval) && pval < 0.05) "**VIOLATION**" else "Pass"
             cat(sprintf("       %-20s: p=%.4f [%s]\n", v, pval, sig), file = result_file, append = TRUE)
        }
    }
    cat("\n", file = result_file, append = TRUE)
  }
  message("Assumption tests complete. Results saved to ", result_file)
}

# ==============================================================================
# 3. PHASE 2: FORWARD AIC LOOP (Linearty Check Only)
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
          r <- cor(1:(length(idx)+1), c(0, coefs[idx]))
          status <- if(r>0.9) "[PASS]" else if(r>0.7) "[WARN]" else "[FAIL]"
          message(sprintf("Linearity Check %-15s: %s (r=%.3f)", v, status, r))
        }
      }
    }
  }
}

run_forward_aic <- function(metadata, data_prefix, target_groups = NULL, output_dir) {
  
  summary_file  <- file.path(output_dir, "aic_results.csv")
  detailed_file <- file.path(output_dir, "clmm_results.csv")
  
  if (!is.null(target_groups)) {
    valid_groups <- intersect(names(metadata), target_groups)
    if (length(valid_groups) == 0) return(list())
    metadata <- metadata[valid_groups]
  }

  flat_metadata <- list()
  for (group in names(metadata)) {
    if (is.list(metadata[[group]]) && !("ind_vars" %in% names(metadata[[group]]))) {
       for (hnum in names(metadata[[group]])) flat_metadata[[hnum]] <- metadata[[group]][[hnum]]
    } else {
       flat_metadata[[group]] <- metadata[[group]]
    }
  }
  metadata <- flat_metadata 
  
  mc <- max(1L, if (is.null(getOption("mc.cores"))) parallel::detectCores(logical=TRUE)-1L else as.integer(getOption("mc.cores")))
  message(sprintf("Using %d parallel worker(s) for AIC loop.", mc))

  for (hnum in names(metadata)) {
    gc() 
    suite_name <- sub("_[0-9]+$", "", hnum)
    suite_dir <- file.path(output_dir, suite_name)
    if (!dir.exists(suite_dir)) dir.create(suite_dir, recursive = TRUE)
    
    log_file <- file.path(suite_dir, paste0(hnum, ".log"))
    log_con <- file(log_file, open = "wt")
    sink(log_con, type = "output", split = TRUE); sink(log_con, type = "message")
    
    tryCatch({
      message("Processing: ", hnum)
      
      ind_vars <- metadata[[hnum]][["ind_vars"]]
      dep_var  <- metadata[[hnum]][["dep_var"]]
      inters   <- metadata[[hnum]][["interactions"]]
      base_predictors <- metadata[[hnum]][["base_predictors"]]
      if (is.null(base_predictors)) base_predictors <- character(0)

      all_vars <- unique(c(ind_vars, base_predictors))
      data_file_path <- file.path(data_prefix, suite_name, paste0(hnum, ".csv"))

      if (!file.exists(data_file_path)) { message("Skipping (No Data)"); sink(type="message"); sink(type="output"); close(log_con); next }
      
      data <- read.csv(data_file_path, check.names = FALSE)
      for (v in all_vars) if (v %in% names(data)) data[[v]] <- cast_by_name(data[[v]], v)
      data[[dep_var]] <- factor(data[[dep_var]], ordered = TRUE)

      do_linear_only <- !grepl("_do_full$", hnum)
      interaction_defs <- build_pairwise_interactions(inters, data, linear_only = do_linear_only)

      build_formula <- function(kept_vars, kept_intrs_strs) {
        main_terms <- vapply(kept_vars, function(v) format_var(v, data, linear_only = do_linear_only), character(1))
        rhs_terms <- c(main_terms, kept_intrs_strs)
        rhs_str   <- paste(rhs_terms, collapse = " + ")
        stats::as.formula(paste(dep_var, "~", if (nzchar(rhs_str)) rhs_str else "1", "+ (1 | topic_condition) + (1 | participant_number) + (1 | data_practice)"))
      }

      kept_vars  <- base_predictors; kept_intrs <- character(0)
      data_fixed <- build_data_sub(data, dep_var, all_vars)
      control <- ordinal::clmm.control(gradTol = 1e-6)
      
      fit <- tryCatch(ordinal::clmm(build_formula(kept_vars, kept_intrs), data = data_fixed, Hess = TRUE, method = "nlminb", control = control), error = function(e) NULL)

      if (is.null(fit)) { message("Base model failed."); sink(type="message"); sink(type="output"); close(log_con); next }
      best_aic <- AIC(fit)
      message(sprintf("Base AIC: %.3f", best_aic))

      improved <- TRUE
      while (improved) {
        improved <- FALSE
        cand_vars <- setdiff(all_vars, kept_vars)
        cand_intrs <- list()
        if (length(interaction_defs) > 0) {
          for (item in interaction_defs) if (!item$term_str %in% kept_intrs && all(item$parents %in% kept_vars)) cand_intrs <- c(cand_intrs, list(item))
        }

        jobs <- list()
        for (v in cand_vars) jobs[[length(jobs)+1]] <- list(type="var", val=v)
        for (i in cand_intrs) jobs[[length(jobs)+1]] <- list(type="intr", val=i$term_str)

        if (length(jobs) == 0) break
        
        trial_results <- parallel::mclapply(jobs, function(job) {
            try_vars <- kept_vars; try_intrs <- kept_intrs
            if (job$type == "var") try_vars <- c(try_vars, job$val) else try_intrs <- c(try_intrs, job$val)
            fit_try <- tryCatch(ordinal::clmm(build_formula(try_vars, try_intrs), data = data_fixed, Hess = TRUE, method = "nlminb", control = control), error = function(e) NULL)
            if (is.null(fit_try)) return(NULL)
            list(job=job, aic=AIC(fit_try), vars=try_vars, intrs=try_intrs)
          }, mc.cores = mc)

        trial_results <- Filter(Negate(is.null), trial_results)
        if (length(trial_results) == 0) break
        
        aics <- vapply(trial_results, function(x) x$aic, numeric(1))
        best_idx <- which.min(aics)
        if (aics[best_idx] < best_aic) {
          best_try <- trial_results[[best_idx]]
          message(sprintf("Add %s: AIC %.3f", if(best_try$job$type=="var") best_try$job$val else best_try$job$val, best_try$aic))
          kept_vars <- best_try$vars; kept_intrs <- best_try$intrs; best_aic <- best_try$aic
          fit <- ordinal::clmm(build_formula(kept_vars, kept_intrs), data = data_fixed, Hess = TRUE, method = "nlminb", control = control)
          improved <- TRUE
        }
      }

      print(summary(fit))
      perform_linearity_check(fit, data_fixed, kept_vars, kept_intrs, dep_var, control)
      
      row_data <- data.frame(
        hypothesis = hnum, AIC = as.numeric(best_aic), nobs = fit$info$nobs, max_grad = fit$info$max.grad,
        k_predictors = length(kept_vars), kept_predictors = paste(kept_vars, collapse = " + "),
        interactions = paste(kept_intrs, collapse = " + "), formula = paste(deparse(build_formula(kept_vars, kept_intrs)), collapse = "")
      )
      write.table(row_data, file = summary_file, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(summary_file))
      
      coefs <- summary(fit)$coefficients
      det_df <- data.frame(Hypothesis=hnum, Term=rownames(coefs), Estimate=coefs[,"Estimate"], Std_Error=coefs[,"Std. Error"], P_Value=coefs[,"Pr(>|z|)"])
      write.table(det_df, file = detailed_file, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(detailed_file))

    }, finally = { sink(type="message"); sink(type="output"); close(log_con) })
  }
}

# ==============================================================================
# 4. DRIVER
# ==============================================================================

# CONFIGURATION FLAGS
DO_ASSUMPTIONS <- TRUE
DO_MODELING    <- TRUE
TARGET_GROUPS  <- c("v3", "v4", "v6")

file_prefix <- "metadata/"
metadata    <- fromJSON(paste0(file_prefix, "metadata.json"))
data_prefix <- "hypothesis_data/"
timestamp   <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
run_dir     <- file.path("results", paste0("run_", timestamp))
if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE)

if (DO_ASSUMPTIONS) {
  message("Starting Assumption Tests...")
  run_all_assumption_tests(metadata, data_prefix, TARGET_GROUPS, run_dir)
}

if (DO_MODELING) {
  message("Starting Forward AIC Modeling...")
  run_forward_aic(metadata, data_prefix, TARGET_GROUPS, run_dir)
}

message("All tasks complete. Results in: ", run_dir)