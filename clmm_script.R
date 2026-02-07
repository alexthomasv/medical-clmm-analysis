# install.packages("ordinal")
# install.packages("insight")
library(jsonlite)
library(ordinal)
library(insight)
library(parallel)

# ---- Formatting/helpers (Unchanged) ----

format_var <- function(var, data, linear_only = TRUE) {
  if (is.ordered(data[[var]]) && isTRUE(linear_only)) {
    sprintf("C(`%s`, contr.poly, 1)", var)   # Linear component only
  } else {
    sprintf("`%s`", var)                      # Full coding
  }
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

  if (is.data.frame(interactions)) {
    interactions <- as.matrix(interactions)
  }

  if (!is.matrix(interactions) || any(dim(interactions) == 0L)) {
    return(character(0))
  }

  vars <- unique(as.vector(interactions))
  vars <- vars[!is.na(vars)]
  vars <- trimws(as.character(vars))
  vars <- vars[nzchar(vars)]
  missing <- setdiff(vars, names(data))
  if (length(missing)) {
    stop("Missing variables in `data`: ", paste(missing, collapse = ", "))
  }

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
        sprintf("C(`%s`, contr.poly, 1)", v)
      } else {
        sprintf("`%s`", v)
      }
    })

    list(
      term_str = paste(term_pieces, collapse = ":"),
      parents  = r
    )
  })

  out_list[!sapply(out_list, is.null)]
}

# ---- Casting helpers (Unchanged) ----

cast_by_name <- function(x, name) {
  lev_usefreq     <- c("Less than once a year","A few times a year",
                       "A few times a month","A few times a week","Daily")
  lev_gender      <- c("Woman","Man","Non-binary", "NO_DATA")
  lev_ethnicity   <- c("White","Asian","Black","Latino/Hispanic","Mixed","Other", "NO_DATA")
  lev_age         <- c("18-20","21-44","45-64","65+", "NO_DATA")
  lev_region      <- c("South","Northeast","West","Midwest","NO_DATA")
  lev_experience_hf <- c("Goodish", "Baddish", "Never Used One")
  lev_experience_med <- c("Goodish", "Baddish", "Never Used One")

  x <- trimws(as.character(x))
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
  bad_nodata <- "NO_DATA"

  if ("gender"         %in% names(df)) df <- df[df$gender         != bad_nodata  & !is.na(df$gender), , drop = FALSE]
  if ("experience_hf"  %in% names(df)) df <- df[df$experience_hf  != bad_never   & !is.na(df$experience_hf), , drop = FALSE]
  if ("experience_med" %in% names(df)) df <- df[df$experience_med != bad_never   & !is.na(df$experience_med), , drop = FALSE]

  df[stats::complete.cases(df[, vars_used, drop = FALSE]), , drop = FALSE]
}

# ---- PARALLEL FORWARD AIC ----

run_forward_aic <- function(metadata, data_prefix, target_groups = NULL, output_dir) {
  
  # Prepare Master Summary Files
  summary_file  <- file.path(output_dir, "aic_results.csv")
  detailed_file <- file.path(output_dir, "clmm_results.csv")
  
  # --- 1. Filter by Top-Level Group ---
  if (!is.null(target_groups)) {
    available_groups <- names(metadata)
    valid_groups <- intersect(available_groups, target_groups)
    if (length(valid_groups) == 0) return(list())
    metadata <- metadata[valid_groups]
  }

  # --- 2. Flatten Metadata Structure ---
  flat_metadata <- list()
  for (group in names(metadata)) {
    if (is.list(metadata[[group]]) && !("ind_vars" %in% names(metadata[[group]]))) {
       for (hnum in names(metadata[[group]])) {
         flat_metadata[[hnum]] <- metadata[[group]][[hnum]]
       }
    } else {
       flat_metadata[[group]] <- metadata[[group]]
    }
  }
  
  metadata <- flat_metadata 
  hnums_to_run <- names(metadata)
  
  message(sprintf("Queued %d hypotheses.\nRun Directory: %s", length(hnums_to_run), output_dir))

  # --- Parallel Setup ---
  results <- list()
  mc_opt      <- getOption("mc.cores")
  mc_detected <- parallel::detectCores(logical = TRUE)
  mc <- max(1L, if (is.null(mc_opt)) mc_detected - 1L else as.integer(mc_opt))
  message(sprintf("Using %d parallel worker(s).", mc))

  for (hnum in hnums_to_run) {
    
    # Force Garbage Collection between hypotheses to clear previous run's memory
    gc() 
    
    suite_name <- sub("_[0-9]+$", "", hnum)
    suite_dir <- file.path(output_dir, suite_name)
    if (!dir.exists(suite_dir)) dir.create(suite_dir, recursive = TRUE)
    
    log_file <- file.path(suite_dir, paste0(hnum, ".log"))
    log_con <- file(log_file, open = "wt")
    sink(log_con, type = "output", split = TRUE)
    sink(log_con, type = "message") 
    
    tryCatch({
      message("\n========================================")
      message(" Processing: ", hnum)
      message("========================================\n")

      ind_vars <- metadata[[hnum]][["ind_vars"]]
      dep_var  <- metadata[[hnum]][["dep_var"]]
      inters   <- metadata[[hnum]][["interactions"]]
      base_predictors <- metadata[[hnum]][["base_predictors"]]
      if (is.null(base_predictors)) base_predictors <- character(0)

      all_vars_to_cast <- unique(c(ind_vars, base_predictors))

      data_file_path <- file.path(data_prefix, suite_name, paste0(hnum, ".csv"))

      if (!file.exists(data_file_path)) {
        message("Skipping ", hnum, " (file not found)")
        sink(type = "message"); sink(type = "output"); close(log_con)
        next
      }
      
      data <- read.csv(data_file_path, check.names = FALSE)

      for (v in all_vars_to_cast) {
        if (v %in% names(data)) data[[v]] <- cast_by_name(data[[v]], v)
      }
      data[[dep_var]] <- factor(data[[dep_var]], ordered = TRUE)

      do_linear_only <- !grepl("_do_full$", hnum)
      interaction_defs <- build_pairwise_interactions(inters, data, linear_only = do_linear_only)

      build_formula <- function(kept_vars, kept_intrs_strs) {
        main_terms <- vapply(kept_vars, function(v) format_var(v, data, linear_only = do_linear_only), character(1))
        rhs_terms <- c(main_terms, kept_intrs_strs)
        rhs_str   <- paste(rhs_terms, collapse = " + ")
        stats::as.formula(paste(dep_var, "~", if (nzchar(rhs_str)) rhs_str else "1", "+ (1 | topic_condition) + (1 | participant_number) + (1 | data_practice)"))
      }

      # ---- STEP 0: BASE MODEL ----
      kept_vars  <- base_predictors
      kept_intrs <- character(0)
      data_fixed <- build_data_sub(data, dep_var, all_vars_to_cast)
      control <- ordinal::clmm.control(gradTol = 1e-6)

      f_curr <- build_formula(kept_vars, kept_intrs)
      print(f_curr) 

      fit <- tryCatch(
        ordinal::clmm(f_curr, data = data_fixed, Hess = TRUE, method = "nlminb", control = control),
        error = function(e) {
          message("Base model failed: ", conditionMessage(e))
          return(NULL)
        }
      )

      if (is.null(fit)) {
         sink(type = "message"); sink(type = "output"); close(log_con)
         next
      }

      best_aic <- AIC(fit)
      message(sprintf("Start (Base) AIC: %.3f", best_aic))

      # ---- FORWARD LOOP ----
      improved <- TRUE
      while (improved) {
        improved <- FALSE
        cand_vars <- setdiff(all_vars_to_cast, kept_vars)
        cand_intrs <- list()
        if (length(interaction_defs) > 0) {
          for (item in interaction_defs) {
            if (item$term_str %in% kept_intrs) next
            if (all(item$parents %in% kept_vars)) cand_intrs <- c(cand_intrs, list(item))
          }
        }

        jobs <- list()
        for (v in cand_vars) jobs[[length(jobs)+1]] <- list(type="var", val=v)
        for (i in cand_intrs) jobs[[length(jobs)+1]] <- list(type="intr", val=i$term_str)

        if (length(jobs) == 0) break

        # Force GC before forking
        gc() 
        
        trial_results <- parallel::mclapply(jobs, function(job) {
            try_vars  <- kept_vars
            try_intrs <- kept_intrs
            if (job$type == "var") try_vars <- c(try_vars, job$val)
            else try_intrs <- c(try_intrs, job$val)

            f_try <- build_formula(try_vars, try_intrs)
            fit_try <- tryCatch(ordinal::clmm(f_try, data = data_fixed, Hess = TRUE, method = "nlminb", control = control), error = function(e) NULL)
            
            if (is.null(fit_try)) return(NULL)
            
            # ### MEMORY FIX: Do NOT return 'fit=fit_try'. 
            # Return only lightweight metadata.
            list(job=job, aic=AIC(fit_try), vars=try_vars, intrs=try_intrs)
            
          }, mc.cores = mc)

        trial_results <- Filter(Negate(is.null), trial_results)
        
        if (length(trial_results) == 0) break
        
        aics <- vapply(trial_results, function(x) as.numeric(x$aic), numeric(1))
        valid_idx <- which(is.finite(aics))

        if (length(valid_idx) == 0) break
        best_try <- trial_results[[valid_idx[which.min(aics[valid_idx])]]]

        if (best_try$aic < best_aic) {
          added_name <- if(best_try$job$type=="var") best_try$job$val else paste0("Intr(", best_try$job$val, ")")
          message(sprintf("Add  %-18s: AIC %.3f -> %.3f", added_name, best_aic, best_try$aic))
          
          kept_vars  <- best_try$vars
          kept_intrs <- best_try$intrs
          best_aic   <- best_try$aic
          
          # ### MEMORY FIX: Re-fit ONLY the winner here.
          # We didn't save the fit object, so we generate it now.
          f_best <- build_formula(kept_vars, kept_intrs)
          fit <- ordinal::clmm(f_best, data = data_fixed, Hess = TRUE, method = "nlminb", control = control)
          
          improved   <- TRUE
        }
      }

      # ---- 1. DUMP FINAL MODEL SUMMARY TO LOG ----
      message("\n--- FINAL MODEL SUMMARY ---")
      print(summary(fit)) 
      message("---------------------------\n")

      # ---- 2. SAVE SUMMARY ROW ----
      final_formula_str <- paste(deparse(build_formula(kept_vars, kept_intrs)), collapse = "")
      row_data <- data.frame(
        hypothesis      = hnum,
        AIC             = as.numeric(best_aic),
        nobs            = as.integer(tryCatch(fit$info$nobs, error = function(e) NA)),
        max_grad        = as.numeric(tryCatch(fit$info$max.grad, error = function(e) NA)),
        k_predictors    = length(kept_vars),
        kept_predictors = paste(kept_vars, collapse = " + "),
        interactions    = if (length(kept_intrs)) paste(kept_intrs, collapse = " + ") else "",
        formula         = final_formula_str,
        stringsAsFactors = FALSE
      )

      write.table(row_data, file = summary_file, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(summary_file))

      # ---- 3. CALCULATE & SAVE DETAILED STATS ----
      tryCatch({
        coef_summary <- summary(fit)$coefficients
        estimates <- coef_summary[, "Estimate"]
        std_errors <- coef_summary[, "Std. Error"]
        
        detailed_df <- data.frame(
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
        
        write.table(detailed_df, file = detailed_file, append = TRUE, sep = ",", row.names = FALSE, col.names = !file.exists(detailed_file))
        message("Saved summary CSVs.")
      }, error = function(e) {
        warning("Failed to save detailed stats for ", hnum, ": ", conditionMessage(e))
      })

      results[[hnum]] <- list(hnum=hnum, aic=best_aic)
      
    }, finally = {
      sink(type = "message")
      sink(type = "output")
      close(log_con)
    })
  }
  invisible(results)
}

# ---- Driver ----
file_prefix <- "metadata/"
file_path   <- paste0(file_prefix, "metadata.json")
metadata    <- fromJSON(file_path)

data_prefix <- "hypothesis_data/"

# 1. Create Timestamped Run Directory
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
run_dir <- file.path("results", paste0("run_", timestamp))

if (!dir.exists(run_dir)) {
  dir.create(run_dir, recursive = TRUE)
}

# 2. Run
res <- run_forward_aic(
    metadata, 
    data_prefix, 
    target_groups = c("v2", "v3", "v4", "v6"),
    output_dir = run_dir
)