source("clmm_script.R")

pass <- 0L
fail <- 0L

assert <- function(expr, msg) {
  result <- tryCatch(eval(expr, envir = parent.frame()), error = function(e) {
    structure(FALSE, error_msg = conditionMessage(e))
  })
  if (isTRUE(result)) {
    pass <<- pass + 1L
    cat(sprintf("  [PASS] %s\n", msg))
  } else {
    fail <<- fail + 1L
    err <- attr(result, "error_msg")
    if (!is.null(err)) {
      cat(sprintf("  [FAIL] %s  (Error: %s)\n", msg, err))
    } else {
      cat(sprintf("  [FAIL] %s  (Got: %s)\n", msg, paste(deparse(result), collapse = "")))
    }
  }
}

assert_error <- function(expr, msg) {
  err <- tryCatch({ eval(expr, envir = parent.frame()); NULL }, error = function(e) e)
  if (!is.null(err)) {
    pass <<- pass + 1L
    cat(sprintf("  [PASS] %s\n", msg))
  } else {
    fail <<- fail + 1L
    cat(sprintf("  [FAIL] %s  (Expected error but succeeded)\n", msg))
  }
}

assert_equal <- function(actual, expected, msg) {
  if (identical(actual, expected)) {
    pass <<- pass + 1L
    cat(sprintf("  [PASS] %s\n", msg))
  } else {
    fail <<- fail + 1L
    cat(sprintf("  [FAIL] %s\n    Expected: %s\n    Got:      %s\n", msg,
                paste(deparse(expected), collapse = ""),
                paste(deparse(actual), collapse = "")))
  }
}

# ===========================================================================
# Synthetic data used across tests
# ===========================================================================
make_test_data <- function() {
  data.frame(
    usefreq_hf     = cast_by_name(c("Daily", "A few times a week", "A few times a month"), "usefreq_hf"),
    experience_hf  = cast_by_name(c("Goodish", "Baddish", "Never Used One"), "experience_hf"),
    age_bracket    = cast_by_name(c("18-20", "21-44", "65+"), "age_bracket"),
    gender         = cast_by_name(c("Woman", "Man", "Non-binary"), "gender"),
    region_broad   = cast_by_name(c("South", "Northeast", "West"), "region_broad"),
    dep            = factor(c("low", "mid", "high"), levels = c("low", "mid", "high"), ordered = TRUE),
    topic_condition     = factor(c("A", "B", "A")),
    participant_number  = factor(c("P1", "P2", "P3")),
    data_practice       = factor(c("dp1", "dp1", "dp2")),
    stringsAsFactors = FALSE
  )
}

# ===========================================================================
cat("\n=== format_var ===\n")
# ===========================================================================

d <- make_test_data()

assert_equal(
  format_var("usefreq_hf", d, linear_only = TRUE, linear_method = "numeric"),
  "scale(as.numeric(`usefreq_hf`))",
  "Ordered + numeric -> scale(as.numeric(...))"
)

assert_equal(
  format_var("usefreq_hf", d, linear_only = TRUE, linear_method = "contrast"),
  "C(`usefreq_hf`, poly, 1)",
  "Ordered + contrast -> C(..., poly, 1)"
)

assert_equal(
  format_var("usefreq_hf", d, linear_only = FALSE),
  "`usefreq_hf`",
  "Ordered + linear_only=FALSE -> backtick-quoted"
)

assert_equal(
  format_var("experience_hf", d, linear_only = TRUE, linear_method = "numeric"),
  "`experience_hf`",
  "Non-ordered -> backtick-quoted (numeric method ignored)"
)

assert_equal(
  format_var("gender", d, linear_only = TRUE, linear_method = "contrast"),
  "`gender`",
  "Non-ordered -> backtick-quoted (contrast method ignored)"
)

# ===========================================================================
cat("\n=== build_interaction_terms ===\n")
# ===========================================================================

assert_equal(
  build_interaction_terms(NULL, d),
  character(0),
  "NULL -> character(0)"
)

assert_equal(
  build_interaction_terms(list(), d),
  character(0),
  "empty list -> character(0)"
)

assert_equal(
  build_interaction_terms("not_a_list", d),
  character(0),
  "non-list -> character(0)"
)

assert_equal(
  build_interaction_terms(list(c("usefreq_hf", "age_bracket")), d, linear_2way = FALSE),
  "`usefreq_hf`:`age_bracket`",
  "2-way non-linear -> backtick-quoted with :"
)

assert_equal(
  build_interaction_terms(list(c("usefreq_hf", "age_bracket")), d, linear_2way = TRUE, linear_method = "numeric"),
  "scale(as.numeric(`usefreq_hf`)):scale(as.numeric(`age_bracket`))",
  "2-way linear numeric -> scale(as.numeric(...)):..."
)

assert_equal(
  build_interaction_terms(list(c("usefreq_hf", "age_bracket", "usefreq_hf")), d, linear_3way = TRUE, linear_method = "numeric"),
  "scale(as.numeric(`usefreq_hf`)):scale(as.numeric(`age_bracket`)):scale(as.numeric(`usefreq_hf`))",
  "3-way linear -> all components linearized"
)

assert_error(
  build_interaction_terms(list(c("usefreq_hf", "gender")), d, linear_2way = TRUE),
  "Non-ordered component in linear interaction -> stop()"
)

# ===========================================================================
cat("\n=== cast_by_name ===\n")
# ===========================================================================

uf <- cast_by_name(c("Daily", "A few times a week"), "usefreq_hf")
assert(is.ordered(uf), "usefreq_hf is ordered")
assert_equal(levels(uf),
  c("Less than once a year", "A few times a year", "A few times a month", "A few times a week", "Daily"),
  "usefreq_hf has 5 levels in correct order"
)

ef <- cast_by_name(c("Goodish", "Baddish"), "experience_hf")
assert(!is.ordered(ef), "experience_hf is NOT ordered")

gf <- cast_by_name(c("Woman", "NO_DATA"), "gender")
assert(is.na(gf[2]), "gender NO_DATA -> NA")

af <- cast_by_name(c("18-20", "NO_DATA"), "age_bracket")
assert(is.na(af[2]), "age_bracket NO_DATA -> NA")

rf <- cast_by_name(c("South", "NO_DATA"), "region_broad")
assert(is.na(rf[2]), "region_broad NO_DATA -> NA")

ab <- cast_by_name(c("18-20", "21-44"), "age_bracket")
assert(is.ordered(ab), "age_bracket is ordered")

gf2 <- cast_by_name(c("Woman", "Man"), "gender")
assert(!is.ordered(gf2), "gender is NOT ordered")

unk <- cast_by_name(c("X", "Y"), "some_unknown_name")
assert(is.ordered(unk), "unknown name -> default ordered factor")

ws <- cast_by_name(c("  Daily  ", " A few times a week "), "usefreq_hf")
assert_equal(as.character(ws), c("Daily", "A few times a week"), "whitespace trimmed")

# ===========================================================================
cat("\n=== build_data_sub ===\n")
# ===========================================================================

big_data <- data.frame(
  dep            = factor(c("low", "mid", "high", "low"), levels = c("low", "mid", "high"), ordered = TRUE),
  usefreq_hf     = cast_by_name(c("Daily", "A few times a week", "A few times a month", "Daily"), "usefreq_hf"),
  age_bracket    = cast_by_name(c("18-20", "21-44", "65+", "18-20"), "age_bracket"),
  extra_col      = c(1, 2, 3, 4),
  topic_condition     = factor(c("A", "B", "A", "B")),
  participant_number  = factor(c("P1", "P2", "P3", "P1")),
  data_practice       = factor(c("dp1", "dp1", "dp2", "dp1")),
  stringsAsFactors = FALSE
)

sub1 <- build_data_sub(big_data, "dep", c("usefreq_hf", "age_bracket"))
assert(!("extra_col" %in% names(sub1)), "build_data_sub drops extra columns")
assert(all(c("dep", "usefreq_hf", "age_bracket", "topic_condition", "participant_number", "data_practice") %in% names(sub1)),
       "build_data_sub keeps dep + needed + control vars")

big_data_na <- big_data
big_data_na$usefreq_hf[2] <- NA
sub_na <- build_data_sub(big_data_na, "dep", c("usefreq_hf", "age_bracket"))
assert_equal(nrow(sub_na), 3L, "build_data_sub drops rows with NA")

# experience_hf filter: the filter looks for "I've never used one" but cast_by_name
# produces "Never Used One" — so no rows are dropped (documents existing behavior)
big_data_exp <- data.frame(
  dep = factor(c("low", "mid"), levels = c("low", "mid", "high"), ordered = TRUE),
  experience_hf = cast_by_name(c("Goodish", "Never Used One"), "experience_hf"),
  topic_condition = factor(c("A", "B")),
  participant_number = factor(c("P1", "P2")),
  data_practice = factor(c("dp1", "dp1")),
  stringsAsFactors = FALSE
)
sub_exp <- build_data_sub(big_data_exp, "dep", c("experience_hf"))
assert_equal(nrow(sub_exp), 2L, "experience_hf filter string mismatch is a no-op on cast data")

# ===========================================================================
cat("\n=== flatten_metadata_wrapper ===\n")
# ===========================================================================

nested <- list(
  group1 = list(
    h1 = list(dep_var = "dv1", base_predictors = list("a")),
    h2 = list(dep_var = "dv2", base_predictors = list("b"))
  ),
  group2 = list(
    h3 = list(dep_var = "dv3", base_predictors = list("c"))
  )
)
flat <- flatten_metadata_wrapper(nested)
assert_equal(sort(names(flat)), c("h1", "h2", "h3"), "nested groups flattened to hypothesis level")

already_flat <- list(
  h1 = list(dep_var = "dv1", base_predictors = list("a"))
)
flat2 <- flatten_metadata_wrapper(already_flat)
assert_equal(names(flat2), "h1", "already-flat input unchanged")

assert_equal(length(flatten_metadata_wrapper(list())), 0L,
             "empty -> empty list")

# ===========================================================================
cat("\n=== get_all_vars_from_hyp ===\n")
# ===========================================================================

hyp1 <- list(
  base_predictors = list("a", "b"),
  candidate_terms_predictors = list("c", "d"),
  base_interactions = list(c("a", "b")),
  candidate_terms_interactions = list(c("c", "d"))
)
vars1 <- get_all_vars_from_hyp(hyp1)
assert(all(c("a", "b", "c", "d") %in% vars1), "all fields populated -> unique union")

hyp2 <- list(
  base_predictors = list("a", "b"),
  candidate_terms_predictors = list("a", "c"),
  base_interactions = list(c("a", "b")),
  candidate_terms_interactions = list(c("a", "c"))
)
vars2 <- get_all_vars_from_hyp(hyp2)
assert_equal(length(vars2), length(unique(vars2)), "duplicates deduplicated")

hyp3 <- list(
  base_predictors = list("x"),
  candidate_terms_predictors = list(),
  base_interactions = list(c("x", "y")),
  candidate_terms_interactions = list(c("x", "y", "z"))
)
vars3 <- get_all_vars_from_hyp(hyp3)
assert(all(c("x", "y", "z") %in% vars3), "interaction vars extracted from nested lists")

# ===========================================================================
cat("\n=== can_add_interaction ===\n")
# ===========================================================================

assert(
  can_add_interaction(c("a", "b"), c("a", "b", "c"), list()),
  "2-way, both mains present -> TRUE"
)

assert(
  !can_add_interaction(c("a", "b"), c("a", "c"), list()),
  "2-way, one main missing -> FALSE"
)

assert(
  can_add_interaction(c("a", "b", "c"), c("a", "b", "c"), list(c("a", "b"), c("a", "c"), c("b", "c"))),
  "3-way, all mains + all 2-ways -> TRUE"
)

assert(
  !can_add_interaction(c("a", "b", "c"), c("a", "b", "c"), list(c("a", "b"), c("a", "c"))),
  "3-way, missing a 2-way sub-interaction -> FALSE"
)

assert(
  !can_add_interaction(c("a", "b", "c"), c("a", "b"), list(c("a", "b"))),
  "3-way, missing a main effect -> FALSE"
)

# Order independence
assert(
  can_add_interaction(c("b", "a", "c"), c("c", "a", "b"), list(c("b", "a"), c("c", "a"), c("c", "b"))),
  "3-way order independence (sorted matching)"
)

# ===========================================================================
cat("\n=== reconstruct_best_models_from_csv ===\n")
# ===========================================================================

tmp_csv <- tempfile(fileext = ".csv")
writeLines(c(
  "hypothesis,AIC,kept_predictors,interactions",
  "h1,100.5,usefreq_hf+age_bracket,",
  "h2,200.3,gender,usefreq_hf:age_bracket"
), tmp_csv)

rec <- reconstruct_best_models_from_csv(tmp_csv)
assert_equal(sort(names(rec)), c("h1", "h2"), "predictors-only CSV -> correct list structure")
assert_equal(rec$h1$kept_vars, c("usefreq_hf", "age_bracket"), "h1 predictors parsed")

assert_equal(rec$h2$kept_intrs_list[[1]], c("usefreq_hf", "age_bracket"),
             "colon-joined interactions parsed to list of vectors")

assert(is.null(reconstruct_best_models_from_csv(tempfile())),
       "nonexistent file -> NULL")

file.remove(tmp_csv)

# ===========================================================================
cat("\n=== get_data_file_path ===\n")
# ===========================================================================

assert_equal(
  get_data_file_path("v6_0", "hypothesis_data/"),
  "hypothesis_data//v6/v6_0.csv",
  "v6_0 -> hypothesis_data//v6/v6_0.csv"
)

assert_equal(
  get_data_file_path("v12_34", "data"),
  "data/v12/v12_34.csv",
  "v12_34 -> data/v12/v12_34.csv"
)

# ===========================================================================
cat("\n=== make_clmm_formula ===\n")
# ===========================================================================

f <- make_clmm_formula("dep", "x + y")
f_str <- deparse(f, width.cutoff = 500)
assert(grepl("dep ~ x \\+ y", f_str), "make_clmm_formula includes dep ~ rhs")
assert(grepl("topic_condition", f_str), "make_clmm_formula includes random effects")
assert(grepl("participant_number", f_str), "make_clmm_formula includes participant_number random effect")
assert(grepl("data_practice", f_str), "make_clmm_formula includes data_practice random effect")

# ===========================================================================
cat("\n=== apply_poly_contrasts ===\n")
# ===========================================================================

d_poly <- data.frame(
  dep = factor(c("low", "mid", "high"), ordered = TRUE),
  ord_var = factor(c("A", "B", "C"), levels = c("A", "B", "C"), ordered = TRUE),
  unord_var = factor(c("X", "Y", "Z")),
  stringsAsFactors = FALSE
)
d_poly <- apply_poly_contrasts(d_poly, "dep")
assert(!is.null(attr(d_poly$ord_var, "contrasts")), "apply_poly_contrasts sets contrasts on ordered vars")
assert(is.null(attr(d_poly$unord_var, "contrasts")), "apply_poly_contrasts skips unordered vars")
assert(is.null(attr(d_poly$dep, "contrasts")), "apply_poly_contrasts skips dep_var")

# ===========================================================================
cat("\n=== load_and_prepare_data ===\n")
# ===========================================================================

tmp_data <- tempfile(fileext = ".csv")
write.csv(data.frame(
  dep = c("low", "mid", "high"),
  usefreq_hf = c("Daily", "A few times a week", "A few times a month"),
  topic_condition = c("A", "B", "A"),
  participant_number = c("P1", "P2", "P3"),
  data_practice = c("dp1", "dp1", "dp2"),
  stringsAsFactors = FALSE
), tmp_data, row.names = FALSE)

res <- load_and_prepare_data(tmp_data, "dep", c("usefreq_hf"))
assert(is.ordered(res$raw$usefreq_hf), "load_and_prepare_data casts vars in raw data")
assert(is.ordered(res$raw$dep), "load_and_prepare_data makes dep_var ordered factor")
assert(is.data.frame(res$fixed), "load_and_prepare_data returns fixed data frame")
file.remove(tmp_data)

# ===========================================================================
# Summary
# ===========================================================================
cat(sprintf("\n%d passed, %d failed\n", pass, fail))
if (fail > 0) quit(status = 1)
