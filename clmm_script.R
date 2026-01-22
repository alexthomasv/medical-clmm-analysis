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
  
  # Normalize to a character matrix WITHOUT recycling shorter rows
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
  
  # Sanity check: all vars exist in data
  vars <- unique(as.vector(interactions))
  vars <- vars[!is.na(vars)]
  vars <- trimws(as.character(vars))
  vars <- vars[nzchar(vars)]
  missing <- setdiff(vars, names(data))
  if (length(missing)) {
    stop("Missing variables in `data`: ", paste(missing, collapse = ", "))
  }
  
  # Ensure factors
  for (v in vars) {
    if (!is.factor(data[[v]])) data[[v]] <- factor(data[[v]], ordered = TRUE)
  }
  
  # Return a list of named descriptors so we can track parents later
  # We return a list: list(term = "string", parents = c("A", "B"))
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
  
  # Filter nulls
  out_list[!sapply(out_list, is.null)]
}

# ---- Your casting helpers (Unchanged) ----

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

run_forward_aic <- function(metadata, data_prefix) {
  results <- list()
  
  # --- how many workers? ---
  mc_opt      <- getOption("mc.cores")
  mc_detected <- parallel::detectCores(logical = TRUE)
  mc <- max(
    1L,
    if (is.null(mc_opt)) mc_detected - 1L else as.integer(mc_opt)
  )
  mc <- 1  # Uncomment to debug sequentially
  message(sprintf("[Forward AIC] Using %d parallel worker(s).", mc))
  
  for (hnum in names(metadata)) {
    message("\n=== Forward AIC for: ", hnum, " ===")
    ind_vars <- metadata[[hnum]][["ind_vars"]]
    dep_var  <- metadata[[hnum]][["dep_var"]]
    inters   <- metadata[[hnum]][["interactions"]]
    base_predictors <- metadata[[hnum]][["base_predictors"]]
    
    # ---- load & cast data ----
    data_file_path <- paste0(data_prefix, hnum, ".csv")
    if (!file.exists(data_file_path)) {
      message("Skipping ", hnum, " (file not found)")
      next
    }
    data <- read.csv(data_file_path, check.names = FALSE)
    for (v in ind_vars) data[[v]] <- cast_by_name(data[[v]], v)
    data[[dep_var]] <- factor(data[[dep_var]], ordered = TRUE)
    
    do_linear_only <- !grepl("_do_full$", hnum)
    
    # Process interactions into a structured list so we can check parents
    # Returns list of objects: { term_str, parents }
    interaction_defs <- build_pairwise_interactions(inters, data, linear_only = do_linear_only)
    
    # ---- formula builder ----
    build_formula <- function(kept_vars, kept_intrs_strs) {
      main_terms <- vapply(
        kept_vars,
        function(v) format_var(v, data, linear_only = do_linear_only),
        character(1)
      )
      
      # Combined RHS
      rhs_terms <- c(main_terms, kept_intrs_strs)
      rhs_str   <- paste(rhs_terms, collapse = " + ")
      
      stats::as.formula(paste(
        dep_var, "~", if (nzchar(rhs_str)) rhs_str else "1",
        "+ (1 | topic_condition) + (1 | participant_number) + (1 | data_practice)"
      ))
    }
    
    # ---- INITIAL STATE: NULL MODEL ----
    kept_vars  <- character(0)
    kept_intrs <- character(0) # Stores the term strings
    
    # We must start with a fit on the minimal data (handling NAs for the full potential set is safer, 
    # but for forward selection usually we subset based on current vars. 
    # To be safe against NAs changing row counts, we can filter based on ALL potential vars first,
    # OR dynamic subsetting. Dynamic subsetting is standard but AIC is only comparable on fixed data.
    # -> We will build data_sub using ALL ind_vars initially to fix the dataset size.
    data_fixed <- build_data_sub(data, dep_var, ind_vars)
    
    control <- ordinal::clmm.control(gradTol = 1e-6)
    
    # Fit Null Model
    f_curr <- build_formula(kept_vars, kept_intrs)
    
    print(f_curr)
    
    fit <- tryCatch(
      ordinal::clmm(f_curr, data = data_fixed, Hess = TRUE, method = "nlminb", control = control),
      error = function(e) {
        warning("Null model failed for ", hnum, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if (is.null(fit)) next
    
    best_aic <- AIC(fit)
    message(sprintf("Start (Null) AIC: %.3f", best_aic))
    
    # ---- FORWARD LOOP ----
    improved <- TRUE
    while (improved) {
      improved <- FALSE
      
      # 1. Identify Candidate Variables (Main Effects)
      #    Only those not yet in kept_vars
      cand_vars <- setdiff(ind_vars, kept_vars)
      
      # 2. Identify Candidate Interactions
      #    Only those not in kept_intrs AND where all parents are in kept_vars (Hierarchy)
      cand_intrs <- list()
      if (length(interaction_defs) > 0) {
        for (item in interaction_defs) {
          # If term already kept, skip
          if (item$term_str %in% kept_intrs) next
          
          # Check hierarchy: are all parents in kept_vars?
          if (all(item$parents %in% kept_vars)) {
            cand_intrs <- c(cand_intrs, list(item))
          }
        }
      }
      
      # Create a unified list of "jobs" for parallel processing
      # Each job is: list(type="var", name="X") or list(type="intr", name="str", obj=...)
      jobs <- list()
      for (v in cand_vars) jobs[[length(jobs)+1]] <- list(type="var", val=v)
      for (i in cand_intrs) jobs[[length(jobs)+1]] <- list(type="intr", val=i$term_str)
      
      if (length(jobs) == 0) {
        message("[Forward AIC] No more candidates available.")
        break
      }
      
      # --- Parallel Candidate Fits ---
      trial_results <- parallel::mclapply(
        jobs,
        function(job) {
          try_kept_vars  <- kept_vars
          try_kept_intrs <- kept_intrs
          
          if (job$type == "var") {
            try_kept_vars <- c(try_kept_vars, job$val)
          } else {
            try_kept_intrs <- c(try_kept_intrs, job$val)
          }
          
          f_try <- build_formula(try_kept_vars, try_kept_intrs)
          
          # Fit
          fit_try <- tryCatch(
            ordinal::clmm(
              f_try, 
              data = data_fixed, # Use fixed data so AIC is comparable
              Hess = TRUE, method = "nlminb", control = control
            ),
            error = function(e) NULL
          )
          
          if (is.null(fit_try)) return(NULL)
          
          list(
            job  = job,
            aic  = AIC(fit_try),
            fit  = fit_try,
            kept_vars = try_kept_vars,
            kept_intrs = try_kept_intrs
          )
        },
        mc.cores = mc
      )
      
      # --- Evaluate Results ---
      trial_results <- Filter(Negate(is.null), trial_results)
      
      # Extract valid AICs
      aics <- vapply(trial_results, function(x) as.numeric(x$aic), numeric(1))
      valid_idx <- which(is.finite(aics))
      
      if (length(valid_idx) == 0) {
        message("[Forward AIC] All candidates failed or produced invalid AIC.")
        break
      }
      
      best_local_idx <- valid_idx[which.min(aics[valid_idx])]
      best_try       <- trial_results[[best_local_idx]]
      
      # Check if this beats the current global best
      if (best_try$aic < best_aic) {
        # Determine log message name
        added_name <- if(best_try$job$type=="var") best_try$job$val else paste0("Intr(", best_try$job$val, ")")
        
        message(sprintf(
          "Add  %-18s: AIC %.3f -> %.3f",
          added_name, best_aic, best_try$aic
        ))
        
        # Update State
        kept_vars  <- best_try$kept_vars
        kept_intrs <- best_try$kept_intrs
        fit        <- best_try$fit
        best_aic   <- best_try$aic
        improved   <- TRUE
      } else {
        message("[Forward AIC] No candidates improved AIC. Best candidate was AIC ", min(aics))
      }
    } # end while
    
    # ---- Finalize ----
    final_formula_str <- paste(deparse(build_formula(kept_vars, kept_intrs)), collapse = "")
    nobs    <- tryCatch(fit$info$nobs,     error = function(e) NA)
    maxgrad <- tryCatch(fit$info$max.grad, error = function(e) NA)
    
    message("Final predictors: ", paste(kept_vars, collapse = ", "))
    if(length(kept_intrs)) message("Final interactions: ", paste(kept_intrs, collapse = " + "))
    message(sprintf("Final AIC: %.3f", best_aic))
    
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
    
  } # end for hnum
  
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
  out[order(out$AIC), , drop = FALSE]
}

# ---- Driver ----
file_prefix <- "blues_r_analysis/metadata/"
file_path   <- paste0(file_prefix, "metadata.json")
metadata    <- fromJSON(file_path)

data_prefix <- "blues_r_analysis/hypothesis_data/"

# Run Forward Selection
res <- run_forward_aic(metadata, data_prefix)
summary_df <- summarize_results(res)
print(summary_df)