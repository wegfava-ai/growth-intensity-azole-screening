# ============================================================
# Final reproducible analysis
# Growth intensity in azole agar screening
#
# PUBLIC REPOSITORY VERSION
# ============================================================
#
# This script is a privacy-preserving reconstruction of the final
# manuscript analysis. It DOES NOT contain or download the original
# study dataset.
#
# REQUIRED LOCAL INPUT:
#   data/analysis_input.csv
#
# See:
#   data/input_template.csv
#   data/input_dictionary.csv
#   data/README.md
#
# IMPORTANT:
# To reproduce the manuscript flow counts, the input should represent
# the full project-level analysis dataset, including records that do not
# have complete paired screening and MIC results. The primary and
# secondary score analyses are performed only on eligible paired records.
#
# Primary analysis:
#   exact growth score (0-3) -> WT/NWT distribution
#   NWT proportion with exact binomial 95% CI
#
# Secondary analyses:
#   performance at score >=1, >=2 and 3
#   citizen/source-cluster bootstrap 95% CI
#   score-based AUC with cluster bootstrap
#   score contrasts with cluster-bootstrap risk differences
#
# Conventional screening-performance metrics are SECONDARY and describe
# only the paired analysis set; they should not be interpreted as
# population-level diagnostic-accuracy estimates.
# ============================================================

# ----------------------------
# 0. SETTINGS
# ----------------------------

INPUT_FILE <- Sys.getenv(
  "LATASP_INPUT_FILE",
  unset = file.path("data", "analysis_input.csv")
)

B <- 5000
SEED <- 20261005

# Optional manuscript count validation.
# Keep FALSE for template/external datasets.
CHECK_PUBLISHED_COUNTS <- FALSE

EXPECTED_COUNTS <- list(
  records = 2779L,
  complete_screening = 2370L,
  incomplete_screening = 409L,
  complete_mic = 370L,
  paired = 238L,
  complete_mic_without_complete_screening = 132L
)

# ----------------------------
# 1. PACKAGES
# ----------------------------

packages <- c(
  "dplyr", "tidyr", "stringr", "readr", "purrr", "tibble",
  "ggplot2", "scales", "pROC", "writexl"
)

missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(packages, library, character.only = TRUE))

options(scipen = 999, dplyr.summarise.inform = FALSE)
set.seed(SEED)

dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("output", "figures"), recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 2. HELPERS
# ----------------------------

parse_lab_number <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x[x %in% c("", "-", "NA", "N/A", "na", "n/a")] <- NA_character_
  x <- stringr::str_replace_all(x, ",", ".")
  suppressWarnings(readr::parse_number(x, locale = readr::locale(decimal_mark = ".")))
}

mic_qualifier <- function(x) {
  x <- stringr::str_trim(as.character(x))
  dplyr::case_when(
    stringr::str_detect(x, "^>=") ~ ">=",
    stringr::str_detect(x, "^>")  ~ ">",
    stringr::str_detect(x, "^<=") ~ "<=",
    stringr::str_detect(x, "^<")  ~ "<",
    TRUE ~ "="
  )
}

class_ecv_1 <- function(raw, value) {
  q <- mic_qualifier(raw)
  dplyr::case_when(
    is.na(value) ~ NA_character_,
    q %in% c(">", ">=") & value >= 1 ~ "NWT",
    q %in% c("<", "<=") & value <= 1 ~ "WT",
    q == "=" & value <= 1 ~ "WT",
    q == "=" & value > 1 ~ "NWT",
    TRUE ~ NA_character_
  )
}

exact_binom_ci <- function(k, n, conf.level = 0.95) {
  if (is.na(n) || n == 0) {
    return(tibble(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  bt <- stats::binom.test(k, n, conf.level = conf.level)
  tibble(
    estimate = k / n,
    lower = unname(bt$conf.int[1]),
    upper = unname(bt$conf.int[2])
  )
}

safe_div <- function(a, b) {
  ifelse(is.na(b) | b == 0, NA_real_, a / b)
}

finite_quantile <- function(x, probs = c(0.025, 0.975)) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  unname(stats::quantile(x, probs = probs, na.rm = TRUE))
}

calc_metrics <- function(truth_nwt, pred_positive) {
  ok <- !is.na(truth_nwt) & !is.na(pred_positive)
  truth_nwt <- as.logical(truth_nwt[ok])
  pred_positive <- as.logical(pred_positive[ok])

  tp <- sum(pred_positive & truth_nwt)
  fp <- sum(pred_positive & !truth_nwt)
  tn <- sum(!pred_positive & !truth_nwt)
  fn <- sum(!pred_positive & truth_nwt)

  sensitivity <- safe_div(tp, tp + fn)
  specificity <- safe_div(tn, tn + fp)
  ppv <- safe_div(tp, tp + fp)
  npv <- safe_div(tn, tn + fn)
  lr_pos <- ifelse(!is.na(sensitivity) & !is.na(specificity) &
                     (1 - specificity) > 0,
                   sensitivity / (1 - specificity), NA_real_)
  lr_neg <- ifelse(!is.na(sensitivity) & !is.na(specificity) &
                     specificity > 0,
                   (1 - sensitivity) / specificity, NA_real_)
  accuracy <- safe_div(tp + tn, tp + fp + tn + fn)
  balanced_accuracy <- ifelse(
    is.na(sensitivity) | is.na(specificity),
    NA_real_,
    (sensitivity + specificity) / 2
  )

  c(
    n = tp + fp + tn + fn,
    tp = tp, fp = fp, tn = tn, fn = fn,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv, npv = npv,
    lr_pos = lr_pos, lr_neg = lr_neg,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy
  )
}

cluster_boot_metrics <- function(data, cutoff, B = B) {
  d <- data |>
    dplyr::filter(!is.na(reference), !is.na(score), !is.na(cluster_id))

  clusters <- unique(d$cluster_id)
  observed <- calc_metrics(d$reference == "NWT", d$score >= cutoff)

  metric_names <- c(
    "sensitivity", "specificity", "ppv", "npv",
    "lr_pos", "lr_neg", "accuracy", "balanced_accuracy"
  )

  boot_mat <- matrix(
    NA_real_, nrow = B, ncol = length(metric_names),
    dimnames = list(NULL, metric_names)
  )

  for (b in seq_len(B)) {
    sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
    idx <- unlist(
      lapply(sampled_clusters, function(cl) which(d$cluster_id == cl)),
      use.names = FALSE
    )
    db <- d[idx, , drop = FALSE]
    m <- calc_metrics(db$reference == "NWT", db$score >= cutoff)
    boot_mat[b, metric_names] <- m[metric_names]
  }

  estimates <- purrr::map_dfr(metric_names, function(metric) {
    ci <- finite_quantile(boot_mat[, metric])
    tibble(
      metric = metric,
      estimate = unname(observed[metric]),
      lower = ci[1],
      upper = ci[2]
    )
  })

  counts <- tibble(
    metric = c("n", "tp", "fp", "tn", "fn"),
    estimate = as.numeric(observed[c("n", "tp", "fp", "tn", "fn")]),
    lower = NA_real_,
    upper = NA_real_
  )

  bind_rows(counts, estimates)
}

calc_auc <- function(data) {
  d <- data |> dplyr::filter(!is.na(reference), !is.na(score))
  if (dplyr::n_distinct(d$reference) < 2 ||
      dplyr::n_distinct(d$score) < 2) return(NA_real_)

  roc_obj <- pROC::roc(
    response = factor(d$reference, levels = c("WT", "NWT")),
    predictor = d$score,
    levels = c("WT", "NWT"),
    direction = "<",
    quiet = TRUE
  )
  as.numeric(pROC::auc(roc_obj))
}

cluster_boot_auc <- function(data, B = B) {
  d <- data |>
    dplyr::filter(!is.na(reference), !is.na(score), !is.na(cluster_id))
  clusters <- unique(d$cluster_id)
  observed <- calc_auc(d)
  boot_auc <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
    idx <- unlist(
      lapply(sampled_clusters, function(cl) which(d$cluster_id == cl)),
      use.names = FALSE
    )
    boot_auc[b] <- calc_auc(d[idx, , drop = FALSE])
  }

  ci <- finite_quantile(boot_auc)
  tibble(estimate = observed, lower = ci[1], upper = ci[2])
}

risk_difference <- function(data, contrast) {
  if (contrast == "Score 1 vs Score 0") {
    d <- data |>
      dplyr::filter(score %in% c(0, 1)) |>
      dplyr::mutate(exposed = score == 1)
  } else if (contrast == "Score >=2 vs Score 0") {
    d <- data |>
      dplyr::filter(score == 0 | score >= 2) |>
      dplyr::mutate(exposed = score >= 2)
  } else {
    stop("Unknown contrast.")
  }

  if (sum(d$exposed, na.rm = TRUE) == 0 ||
      sum(!d$exposed, na.rm = TRUE) == 0) {
    return(c(
      p_exposed = NA_real_,
      p_score0 = NA_real_,
      risk_difference = NA_real_
    ))
  }

  p_exposed <- mean(d$reference[d$exposed] == "NWT", na.rm = TRUE)
  p_score0 <- mean(d$reference[!d$exposed] == "NWT", na.rm = TRUE)

  c(
    p_exposed = p_exposed,
    p_score0 = p_score0,
    risk_difference = p_exposed - p_score0
  )
}

cluster_boot_rd <- function(data, contrast, B = B) {
  d <- data |>
    dplyr::filter(!is.na(reference), !is.na(score), !is.na(cluster_id))
  clusters <- unique(d$cluster_id)
  observed <- risk_difference(d, contrast)
  boot_rd <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
    idx <- unlist(
      lapply(sampled_clusters, function(cl) which(d$cluster_id == cl)),
      use.names = FALSE
    )
    tmp <- risk_difference(d[idx, , drop = FALSE], contrast)
    boot_rd[b] <- tmp["risk_difference"]
  }

  ci <- finite_quantile(boot_rd)
  tibble(
    contrast = contrast,
    p_exposed = unname(observed["p_exposed"]),
    p_score0 = unname(observed["p_score0"]),
    risk_difference = unname(observed["risk_difference"]),
    lower = ci[1],
    upper = ci[2]
  )
}

check_count <- function(observed, expected, label) {
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(
      "Published-count validation failed for ", label,
      ". Observed=", observed, "; expected=", expected, "."
    )
  }
}

# ----------------------------
# 3. IMPORT PUBLIC ANALYSIS INPUT
# ----------------------------

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input file not found: ", INPUT_FILE, "\n",
    "Create data/analysis_input.csv using data/input_template.csv and ",
    "data/input_dictionary.csv. The real study data are intentionally ",
    "not included in this repository."
  )
}

raw <- readr::read_csv(
  INPUT_FILE,
  col_types = readr::cols(.default = readr::col_character()),
  na = c("", "NA", "N/A")
)

required_cols <- c(
  "isolate_id", "cluster_id", "country",
  "score_itc", "score_vrc", "score_pos",
  "mic_itc", "mic_vrc", "mic_pos"
)

missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Missing required input columns: ", paste(missing_cols, collapse = ", "))
}

if (anyDuplicated(raw$isolate_id)) {
  stop("isolate_id must be unique. Duplicate isolate IDs were detected.")
}

dat <- raw |>
  dplyr::mutate(
    isolate_id = stringr::str_squish(isolate_id),
    cluster_id = stringr::str_squish(cluster_id),
    country = stringr::str_squish(country),

    score_itc_raw = score_itc,
    score_vrc_raw = score_vrc,
    score_pos_raw = score_pos,
    mic_itc_raw = mic_itc,
    mic_vrc_raw = mic_vrc,
    mic_pos_raw = mic_pos,

    score_itc = parse_lab_number(score_itc),
    score_vrc = parse_lab_number(score_vrc),
    score_pos = parse_lab_number(score_pos),
    mic_itc = parse_lab_number(mic_itc),
    mic_vrc = parse_lab_number(mic_vrc),
    mic_pos = parse_lab_number(mic_pos)
  )

bad_scores <- dat |>
  dplyr::filter(
    (!is.na(score_itc) & !score_itc %in% 0:3) |
    (!is.na(score_vrc) & !score_vrc %in% 0:3) |
    (!is.na(score_pos) & !score_pos %in% 0:3)
  )

if (nrow(bad_scores) > 0) {
  stop("Screening scores outside the allowed 0-3 range were detected.")
}

# ----------------------------
# 4. REFERENCE CLASSIFICATION
# ----------------------------

dat <- dat |>
  dplyr::mutate(
    itc_class = class_ecv_1(mic_itc_raw, mic_itc),
    vrc_class = class_ecv_1(mic_vrc_raw, mic_vrc),
    pos_qualifier = mic_qualifier(mic_pos_raw),
    pos_class = dplyr::case_when(
      is.na(mic_pos) ~ NA_character_,
      pos_qualifier %in% c(">", ">=") & mic_pos >= 0.25 ~ "NWT",
      pos_qualifier %in% c("<", "<=") & mic_pos <= 0.125 ~ "WT",
      pos_qualifier == "=" & mic_pos <= 0.125 ~ "WT",
      pos_qualifier == "=" & mic_pos > 0.25 ~ "NWT",
      pos_qualifier == "=" & mic_pos == 0.25 & itc_class == "WT" ~ "WT",
      pos_qualifier == "=" & mic_pos == 0.25 & itc_class == "NWT" ~ "NWT",
      TRUE ~ NA_character_
    ),
    screen_complete3 = !is.na(score_itc) & !is.na(score_vrc) & !is.na(score_pos),
    mic_complete3 = !is.na(mic_itc) & !is.na(mic_vrc) & !is.na(mic_pos),
    refclass_complete3 = !is.na(itc_class) & !is.na(vrc_class) & !is.na(pos_class)
  )

# ----------------------------
# 5. STUDY FLOW COUNTS
# ----------------------------

# These aggregate counts reproduce the manuscript study-flow summary.
flow_counts <- tibble(
  step = c(
    "Project records with isolate ID",
    "Complete 3-azole screening",
    "No complete 3-azole screening",
    "Complete 3-azole MIC",
    "Complete 3-azole MIC without complete 3-azole screening",
    "Final paired complete screening + MIC"
  ),
  n = c(
    nrow(dat),
    sum(dat$screen_complete3),
    sum(!dat$screen_complete3),
    sum(dat$mic_complete3),
    sum(dat$mic_complete3 & !dat$screen_complete3),
    sum(dat$screen_complete3 & dat$mic_complete3 & dat$refclass_complete3)
  )
)

readr::write_csv(flow_counts, "output/tables/study_flow_counts.csv")

if (CHECK_PUBLISHED_COUNTS) {
  check_count(nrow(dat), EXPECTED_COUNTS$records, "records")
  check_count(sum(dat$screen_complete3), EXPECTED_COUNTS$complete_screening, "complete screening")
  check_count(sum(!dat$screen_complete3), EXPECTED_COUNTS$incomplete_screening, "incomplete screening")
  check_count(sum(dat$mic_complete3), EXPECTED_COUNTS$complete_mic, "complete MIC")
  check_count(
    sum(dat$mic_complete3 & !dat$screen_complete3),
    EXPECTED_COUNTS$complete_mic_without_complete_screening,
    "complete MIC without complete screening"
  )
  check_count(
    sum(dat$screen_complete3 & dat$mic_complete3 & dat$refclass_complete3),
    EXPECTED_COUNTS$paired,
    "paired"
  )
}

# ----------------------------
# 6. PAIRED ANALYSIS DATASET
# ----------------------------

paired <- dat |>
  dplyr::filter(screen_complete3, mic_complete3, refclass_complete3) |>
  dplyr::mutate(
    overall_class = dplyr::case_when(
      itc_class == "NWT" | vrc_class == "NWT" | pos_class == "NWT" ~ "NWT",
      itc_class == "WT" & vrc_class == "WT" & pos_class == "WT" ~ "WT",
      TRUE ~ NA_character_
    )
  )

paired_long <- dplyr::bind_rows(
  paired |> dplyr::transmute(
    isolate_id, cluster_id, country,
    azole = "Itraconazole", score = as.integer(score_itc),
    mic = mic_itc, mic_raw = mic_itc_raw, reference = itc_class
  ),
  paired |> dplyr::transmute(
    isolate_id, cluster_id, country,
    azole = "Voriconazole", score = as.integer(score_vrc),
    mic = mic_vrc, mic_raw = mic_vrc_raw, reference = vrc_class
  ),
  paired |> dplyr::transmute(
    isolate_id, cluster_id, country,
    azole = "Posaconazole", score = as.integer(score_pos),
    mic = mic_pos, mic_raw = mic_pos_raw, reference = pos_class
  )
) |>
  dplyr::mutate(
    azole = factor(
      azole,
      levels = c("Itraconazole", "Voriconazole", "Posaconazole")
    ),
    reference = factor(reference, levels = c("WT", "NWT"))
  )

# ----------------------------
# 7. COUNTRY TABLE
# ----------------------------

table1_country <- paired |>
  dplyr::count(country, overall_class, name = "n") |>
  tidyr::complete(
    country,
    overall_class = c("WT", "NWT"),
    fill = list(n = 0)
  ) |>
  tidyr::pivot_wider(
    names_from = overall_class,
    values_from = n,
    values_fill = 0
  ) |>
  dplyr::mutate(
    Total = WT + NWT,
    NWT_percent = 100 * NWT / Total
  ) |>
  dplyr::arrange(dplyr::desc(Total))

readr::write_csv(
  table1_country,
  "output/tables/Table1_country_global_WT_NWT.csv"
)

# ----------------------------
# 8. PRIMARY EXACT-SCORE ANALYSIS
# ----------------------------

score_reference_counts <- paired_long |>
  dplyr::count(azole, score, reference, name = "n") |>
  tidyr::complete(
    azole,
    score = 0:3,
    reference = factor(c("WT", "NWT"), levels = c("WT", "NWT")),
    fill = list(n = 0)
  ) |>
  dplyr::arrange(azole, score, reference)

readr::write_csv(
  score_reference_counts,
  "output/tables/Supplement_score_x_WT_NWT_counts.csv"
)

nwt_by_score <- paired_long |>
  dplyr::group_by(azole, score) |>
  dplyr::summarise(
    n = dplyr::n(),
    n_nwt = sum(reference == "NWT"),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    ci = purrr::map2(n_nwt, n, exact_binom_ci)
  ) |>
  tidyr::unnest(ci) |>
  dplyr::rename(
    proportion_nwt = estimate,
    ci_lower = lower,
    ci_upper = upper
  ) |>
  dplyr::arrange(azole, score)

readr::write_csv(
  nwt_by_score,
  "output/tables/Primary_NWT_proportion_by_score_exact95CI.csv"
)

# ----------------------------
# 9. SECONDARY PERFORMANCE
# ----------------------------

cutoffs <- c(1, 2, 3)

performance_cluster <- paired_long |>
  split(.$azole) |>
  purrr::imap_dfr(function(d_azole, azole_name) {
    purrr::map_dfr(cutoffs, function(cutoff) {
      cluster_boot_metrics(d_azole, cutoff = cutoff, B = B) |>
        dplyr::mutate(azole = azole_name, cutoff = cutoff, .before = 1)
    })
  })

readr::write_csv(
  performance_cluster,
  "output/tables/Secondary_performance_cluster_bootstrap.csv"
)

performance_wide <- performance_cluster |>
  dplyr::filter(!metric %in% c("n", "tp", "fp", "tn", "fn")) |>
  dplyr::mutate(
    formatted = dplyr::case_when(
      metric %in% c(
        "sensitivity", "specificity", "ppv", "npv",
        "accuracy", "balanced_accuracy"
      ) ~ sprintf(
        "%.1f%% (%.1f-%.1f)",
        100 * estimate, 100 * lower, 100 * upper
      ),
      TRUE ~ sprintf("%.2f (%.2f-%.2f)", estimate, lower, upper)
    )
  ) |>
  dplyr::select(azole, cutoff, metric, formatted) |>
  tidyr::pivot_wider(names_from = metric, values_from = formatted)

readr::write_csv(
  performance_wide,
  "output/tables/Secondary_performance_cluster_bootstrap_formatted.csv"
)

auc_cluster <- paired_long |>
  split(.$azole) |>
  purrr::imap_dfr(function(d_azole, azole_name) {
    cluster_boot_auc(d_azole, B = B) |>
      dplyr::mutate(azole = azole_name, .before = 1)
  })

readr::write_csv(
  auc_cluster,
  "output/tables/Secondary_AUC_cluster_bootstrap.csv"
)

# ----------------------------
# 10. SCORE CONTRASTS
# ----------------------------

contrasts <- c("Score 1 vs Score 0", "Score >=2 vs Score 0")

score_contrasts <- paired_long |>
  split(.$azole) |>
  purrr::imap_dfr(function(d_azole, azole_name) {
    purrr::map_dfr(contrasts, function(contrast_name) {
      cluster_boot_rd(d_azole, contrast_name, B = B) |>
        dplyr::mutate(azole = azole_name, .before = 1)
    })
  })

readr::write_csv(
  score_contrasts,
  "output/tables/Score_contrasts_cluster_bootstrap_RD.csv"
)

# ----------------------------
# 11. DISCORDANCE ANALYSIS
# ----------------------------

discordance_all <- paired_long |>
  dplyr::mutate(
    positive_ge1 = score >= 1,
    truth_nwt = reference == "NWT",
    classification_ge1 = dplyr::case_when(
      positive_ge1 & truth_nwt ~ "TP",
      positive_ge1 & !truth_nwt ~ "FP",
      !positive_ge1 & !truth_nwt ~ "TN",
      !positive_ge1 & truth_nwt ~ "FN",
      TRUE ~ NA_character_
    )
  )

discordance_summary <- discordance_all |>
  dplyr::count(azole, classification_ge1, name = "n") |>
  tidyr::complete(
    azole,
    classification_ge1 = c("TP", "FP", "TN", "FN"),
    fill = list(n = 0)
  )

readr::write_csv(
  discordance_summary,
  "output/tables/Discordance_summary_cutoff_ge1.csv"
)

score1_fp <- discordance_all |>
  dplyr::group_by(azole) |>
  dplyr::summarise(
    total_false_positives_ge1 = sum(classification_ge1 == "FP"),
    score1_false_positives = sum(score == 1 & reference == "WT"),
    proportion_of_false_positives_from_score1 =
      safe_div(score1_false_positives, total_false_positives_ge1),
    .groups = "drop"
  )

readr::write_csv(
  score1_fp,
  "output/tables/Score1_contribution_to_false_positives.csv"
)

# Public repository version deliberately does NOT export isolate-level
# discordance rows, to avoid encouraging accidental release of record-level
# study data. Aggregate tables above fully reproduce manuscript summaries.

# ----------------------------
# 12. AUDITS AND CONSOLIDATED WORKBOOK
# ----------------------------

reference_counts <- tibble(
  azole = c("Itraconazole", "Voriconazole", "Posaconazole"),
  WT = c(
    sum(paired$itc_class == "WT"),
    sum(paired$vrc_class == "WT"),
    sum(paired$pos_class == "WT")
  ),
  NWT = c(
    sum(paired$itc_class == "NWT"),
    sum(paired$vrc_class == "NWT"),
    sum(paired$pos_class == "NWT")
  )
) |>
  dplyr::mutate(Total = WT + NWT)

score_counts <- dplyr::bind_rows(
  paired |> dplyr::count(score = score_itc) |>
    dplyr::mutate(azole = "Itraconazole", .before = 1),
  paired |> dplyr::count(score = score_vrc) |>
    dplyr::mutate(azole = "Voriconazole", .before = 1),
  paired |> dplyr::count(score = score_pos) |>
    dplyr::mutate(azole = "Posaconazole", .before = 1)
) |>
  dplyr::arrange(azole, score)

readr::write_csv(reference_counts, "output/tables/reference_classification_counts.csv")
readr::write_csv(score_counts, "output/tables/exact_score_counts.csv")

final_tables <- list(
  study_flow = flow_counts,
  Table1_country = table1_country,
  score_x_reference = score_reference_counts,
  NWT_by_score_exactCI = nwt_by_score,
  score_contrasts_RD = score_contrasts,
  score1_FP_contribution = score1_fp,
  perf_cluster_boot = performance_cluster,
  perf_formatted = performance_wide,
  cross_classification_ge1 = discordance_summary,
  AUC_cluster_boot = auc_cluster,
  reference_counts = reference_counts,
  score_counts = score_counts
)

writexl::write_xlsx(
  final_tables,
  path = "output/LATAsp_reproducibility_tables.xlsx"
)

log_lines <- c(
  "Public reproducibility analysis",
  paste0("Run date: ", Sys.Date()),
  paste0("Input file: ", INPUT_FILE),
  paste0("Bootstrap replicates: ", B),
  paste0("Bootstrap cluster: cluster_id"),
  paste0("Total records: ", nrow(dat)),
  paste0("Complete 3-azole screening: ", sum(dat$screen_complete3)),
  paste0("Complete 3-azole MIC: ", sum(dat$mic_complete3)),
  paste0("Final paired dataset: ", nrow(paired)),
  "",
  "Primary: exact score 0-3 x WT/NWT with exact binomial 95% CI.",
  "Secondary: screening performance, AUC and score contrasts with cluster bootstrap.",
  "Secondary performance estimates describe the paired analysis set only."
)

writeLines(log_lines, "output/analysis_log.txt")
capture.output(sessionInfo(), file = "output/sessionInfo.txt")

message("Final analysis completed. Outputs are in output/.")
