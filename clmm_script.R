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

# UPDATED: Handles subdirectories for data files
run_forward_aic <- function(metadata, data_prefix, target_groups = NULL, output_file = "aic_results_incremental.csv") {
  
  # --- 1. Filter by Top-Level Group ---
  if (!is.null(target_groups)) {
    available_groups <- names(metadata)
    missing <- setdiff(target_groups, available_groups)
    if (length(missing) > 0) {
      warning("The following groups were not found in metadata: ", paste(missing, collapse = ", "))
    }
    
    valid_groups <- intersect(available_groups, target_groups)
    if (length(valid_groups) == 0) {
      message("No valid groups selected. Exiting.")
      return(list())
    }
    
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
  
  message(sprintf("[Forward AIC] Selected Groups: %s", if(is.null(target_groups)) "ALL" else paste(target_groups, collapse=", ")))
  message(sprintf("[Forward AIC] Queued %d hypotheses: %s", length(hnums_to_run), paste(hnums_to_run, collapse=", ")))
  message(sprintf("[Forward AIC] Appending results to: %s", output_file))

  # --- Parallel Setup ---
  results <- list()
  mc_opt      <- getOption("mc.cores")
  mc_detected <- parallel::detectCores(logical = TRUE)
  mc <- max(1L, if (is.null(mc_opt)) mc_detected - 1L else as.integer(mc_opt))
  message(sprintf("[Forward AIC] Using %d parallel worker(s).", mc))

  for (hnum in hnums_to_run) {
    message("\n=== Forward AIC for: ", hnum, " ===")

    ind_vars <- metadata[[hnum]][["ind_vars"]]
    dep_var  <- metadata[[hnum]][["dep_var"]]
    inters   <- metadata[[hnum]][["interactions"]]

    # 1. SETUP: Base Predictors
    base_predictors <- metadata[[hnum]][["base_predictors"]]
    if (is.null(base_predictors)) base_predictors <- character(0)

    all_vars_to_cast <- unique(c(ind_vars, base_predictors))

    # ---- load & cast data ----
    
    # NEW LOGIC: Calculate subdirectory based on hnum
    # This assumes hnum is like "v1_0" (parent "v1") or "v1_dummy_0" (parent "v1_dummy")
    # We strip the last underscore and the digits following it to get the parent folder.
    parent_folder <- sub("_[0-9]+$", "", hnum)
    
    # Construct path: hypothesis_data/v1/v1_0.csv
    data_file_path <- file.path(data_prefix, parent_folder, paste0(hnum, ".csv"))

    if (!file.exists(data_file_path)) {
      message("Skipping ", hnum, " (file not found at: ", data_file_path, ")")
      next
    }
    data <- read.csv(data_file_path, check.names = FALSE)

    for (v in all_vars_to_cast) {
      if (v %in% names(data)) {
        data[[v]] <- cast_by_name(data[[v]], v)
      } else {
        warning(sprintf("Variable '%s' listed in metadata but missing from CSV.", v))
      }
    }
    data[[dep_var]] <- factor(data[[dep_var]], ordered = TRUE)

    do_linear_only <- !grepl("_do_full$", hnum)

    # Process interactions
    interaction_defs <- build_pairwise_interactions(inters, data, linear_only = do_linear_only)

    # ---- formula builder ----
    build_formula <- function(kept_vars, kept_intrs_strs) {
      main_terms <- vapply(
        kept_vars,
        function(v) format_var(v, data, linear_only = do_linear_only),
        character(1)
      )

      rhs_terms <- c(main_terms, kept_intrs_strs)
      rhs_str   <- paste(rhs_terms, collapse = " + ")

      stats::as.formula(paste(
        dep_var, "~", if (nzchar(rhs_str)) rhs_str else "1",
        "+ (1 | topic_condition) + (1 | participant_number) + (1 | data_practice)"
      ))
    }

    # ---- STEP 0: START DIRECTLY WITH BASE MODEL ----
    kept_vars  <- base_predictors
    kept_intrs <- character(0)
    data_fixed <- build_data_sub(data, dep_var, all_vars_to_cast)
    control <- ordinal::clmm.control(gradTol = 1e-6)

    f_curr <- build_formula(kept_vars, kept_intrs)
    print(f_curr)

    fit <- tryCatch(
      ordinal::clmm(f_curr, data = data_fixed, Hess = TRUE, method = "nlminb", control = control),
      error = function(e) {
        warning("Base model failed for ", hnum, ": ", conditionMessage(e))
        return(NULL)
      }
    )

    if (is.null(fit)) next

    best_aic <- AIC(fit)
    message(sprintf("Start (Base) AIC: %.3f | Predictors: %s", best_aic, paste(kept_vars, collapse=", ")))

    # ---- FORWARD LOOP ----
    improved <- TRUE
    while (improved) {
      improved <- FALSE

      cand_vars <- setdiff(all_vars_to_cast, kept_vars)
      cand_intrs <- list()
      if (length(interaction_defs) > 0) {
        for (item in interaction_defs) {
          if (item$term_str %in% kept_intrs) next
          if (all(item$parents %in% kept_vars)) {
            cand_intrs <- c(cand_intrs, list(item))
          }
        }
      }

      jobs <- list()
      for (v in cand_vars) jobs[[length(jobs)+1]] <- list(type="var", val=v)
      for (i in cand_intrs) jobs[[length(jobs)+1]] <- list(type="intr", val=i$term_str)

      if (length(jobs) == 0) break

      trial_results <- parallel::mclapply(
        jobs,
        function(job) {
          try_kept_vars  <- kept_vars
          try_kept_intrs <- kept_intrs
          if (job$type == "var") try_kept_vars <- c(try_kept_vars, job$val)
          else try_kept_intrs <- c(try_kept_intrs, job$val)

          f_try <- build_formula(try_kept_vars, try_kept_intrs)
          fit_try <- tryCatch(
            ordinal::clmm(f_try, data = data_fixed, Hess = TRUE, method = "nlminb", control = control),
            error = function(e) NULL
          )
          if (is.null(fit_try)) return(NULL)
          list(job=job, aic=AIC(fit_try), fit=fit_try, kept_vars=try_kept_vars, kept_intrs=try_kept_intrs)
        },
        mc.cores = mc
      )

      trial_results <- Filter(Negate(is.null), trial_results)
      aics <- vapply(trial_results, function(x) as.numeric(x$aic), numeric(1))
      valid_idx <- which(is.finite(aics))

      if (length(valid_idx) == 0) break

      best_local_idx <- valid_idx[which.min(aics[valid_idx])]
      best_try       <- trial_results[[best_local_idx]]

      if (best_try$aic < best_aic) {
        added_name <- if(best_try$job$type=="var") best_try$job$val else paste0("Intr(", best_try$job$val, ")")
        message(sprintf("Add  %-18s: AIC %.3f -> %.3f", added_name, best_aic, best_try$aic))
        kept_vars  <- best_try$kept_vars
        kept_intrs <- best_try$kept_intrs
        fit        <- best_try$fit
        best_aic   <- best_try$aic
        improved   <- TRUE
      }
    }

    # ---- Final Output Calculation ----
    final_formula_str <- paste(deparse(build_formula(kept_vars, kept_intrs)), collapse = "")
    nobs    <- tryCatch(fit$info$nobs,     error = function(e) NA)
    maxgrad <- tryCatch(fit$info$max.grad, error = function(e) NA)

    # --- SAVE INCREMENTALLY TO FILE ---
    row_data <- data.frame(
      hypothesis      = hnum,
      AIC             = as.numeric(best_aic),
      nobs            = as.integer(nobs),
      max_grad        = as.numeric(maxgrad),
      k_predictors    = length(kept_vars),
      kept_predictors = paste(kept_vars, collapse = " + "),
      interactions    = if (length(kept_intrs)) paste(kept_intrs, collapse = " + ") else "",
      formula         = final_formula_str,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    write.table(
      row_data, 
      file = output_file, 
      append = file.exists(output_file), 
      sep = ",", 
      row.names = FALSE, 
      col.names = !file.exists(output_file)
    )
    
    results[[hnum]] <- list(
      hnum            = hnum,
      final_fit       = fit,
      kept_predictors = kept_vars,
      kept_interactions = kept_intrs,
      AIC             = as.numeric(best_aic),
      nobs            = as.integer(nobs),
      max_grad        = as.numeric(maxgrad),
      formula         = final_formula_str
    )
  }

  invisible(results)
}

summarize_results <- function(results_list) {
  rows <- lapply(results_list, function(x) {
    data.frame(
      hypothesis      = x$hnum,
      AIC             = x$AIC,
      nobs            = x$nobs,
      max_grad        = x$max_grad,
      k_predictors    = length(x$kept_predictors),
      kept_predictors = paste(x$kept_predictors, collapse = " + "),
      interactions    = if (length(x$kept_interactions)) paste(x$kept_interactions, collapse = " + ") else "",
      formula         = x$formula,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out[order(out$AIC), , drop = FALSE]
}

# ---- Driver ----
file_prefix <- "metadata/"
file_path   <- paste0(file_prefix, "metadata.json")
metadata    <- fromJSON(file_path)

data_prefix <- "hypothesis_data/"
output_csv <- "aic_results.csv"

# Optional: Clear old results
if (file.exists(output_csv)) file.remove(output_csv)

res <- run_forward_aic(
    metadata, 
    data_prefix, 
    target_groups = c("v1", "v1_dummy", "v2", "v3", "v4", "v6"),
    output_file = output_csv
)

summary_df <- summarize_results(res)
print(summary_df)