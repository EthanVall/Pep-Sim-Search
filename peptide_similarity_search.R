# ─────────────────────────────────────────────────────────────────────────────
# Peptide Similarity Search — Interactive Version
# ─────────────────────────────────────────────────────────────────────────────
# Usage:
#   source("peptide_similarity_search.R")
#   results <- run_peptide_search(LMD)
#
# The function will prompt you to enter:
#   1. The peptide sequence column name
#   2. The normal/reference intensity column name
#   3. Your query peptide sequence
#   4. How many top results you want
# ─────────────────────────────────────────────────────────────────────────────

# ── 1. Amino acid physicochemical groups ─────────────────────────────────────
AA_GROUPS <- list(
  "1_small_nonpolar"  = c("G", "A", "V", "L", "I", "P"),
  "2_aromatic"        = c("F", "Y", "W"),
  "3_polar_uncharged" = c("S", "T", "C", "M", "N", "Q"),
  "4_positive"        = c("R", "K", "H"),
  "5_negative"        = c("D", "E")
)

AA_TO_GROUP <- setNames(
  unlist(lapply(seq_along(AA_GROUPS), function(i) rep(i, length(AA_GROUPS[[i]])))),
  unlist(AA_GROUPS)
)

# ── 2. Per-position scoring ───────────────────────────────────────────────────
score_position <- function(aa_i, aa_j) {
  if (aa_i == aa_j) return(2L)
  gi <- AA_TO_GROUP[aa_i]
  gj <- AA_TO_GROUP[aa_j]
  if (!is.na(gi) && !is.na(gj) && gi == gj) return(1L)
  return(0L)
}

# ── 3. Score one candidate against query ─────────────────────────────────────
compute_scores <- function(query_aas, candidate) {
  cand_aas <- strsplit(candidate, "")[[1]]
  q_len    <- length(query_aas)
  c_len    <- length(cand_aas)
  overlap  <- min(q_len, c_len)

  pos_scores   <- mapply(score_position, query_aas[1:overlap], cand_aas[1:overlap])
  raw_sim      <- sum(pos_scores)
  len_diff     <- abs(q_len - c_len)
  penalty      <- len_diff * 0.5
  max_possible <- q_len * 2L

  norm_sim    <- round(max(0, (raw_sim  - penalty) / max_possible), 4)
  n_identical <- sum(pos_scores == 2L)
  norm_id     <- round(max(0, (n_identical - penalty) / q_len), 4)

  list(
    norm_similarity = norm_sim,
    raw_similarity  = raw_sim,
    norm_identity   = norm_id,
    n_identical     = n_identical,
    len_penalty     = penalty,
    len_diff        = len_diff,
    cand_length     = c_len,
    position_scores = paste(pos_scores, collapse = "-")
  )
}

# ── 4. Helper: prompt user with validation ────────────────────────────────────
prompt_colname <- function(prompt_text, valid_cols) {
  repeat {
    cat("\nAvailable columns:\n")
    # Print in rows of 4 for readability
    chunks <- split(valid_cols, ceiling(seq_along(valid_cols) / 4))
    for (chunk in chunks) {
      cat(sprintf("  %-35s", chunk), "\n", sep = "")
    }
    cat("\n", prompt_text, ": ", sep = "")
    input <- trimws(readLines(con = stdin(), n = 1))
    if (input %in% valid_cols) return(input)
    cat(sprintf('\n  [!] "%s" not found. Please enter an exact column name from the list above.\n', input))
  }
}

prompt_text_input <- function(prompt_text, validator = NULL, error_msg = "Invalid input.") {
  repeat {
    cat(prompt_text, ": ", sep = "")
    input <- trimws(readLines(con = stdin(), n = 1))
    if (is.null(validator) || validator(input)) return(input)
    cat(paste0("\n  [!] ", error_msg, "\n\n"))
  }
}

# ── 5. Main interactive function ──────────────────────────────────────────────
run_peptide_search <- function(df) {

  stopifnot(is.data.frame(df))

  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║          Peptide Similarity Search — Setup               ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")

  all_cols <- colnames(df)

  # ── Prompt 1: peptide sequence column ──────────────────────────────────────
  cat("\n[1/4] Which column contains the PEPTIDE SEQUENCES?\n")
  pep_col <- prompt_colname("Enter peptide sequence column name", all_cols)
  cat(sprintf('  -> Using "%s" as peptide sequence column.\n', pep_col))

  # ── Prompt 2: normal/reference intensity column ────────────────────────────
  cat("\n[2/4] Which column contains the NORMAL/REFERENCE INTENSITY values?\n")
  cat("      (Only peptides with intensity > 0 in this column will be searched)\n")
  norm_col <- prompt_colname("Enter normal intensity column name", all_cols)
  cat(sprintf('  -> Using "%s" as the reference intensity column.\n', norm_col))

  # ── Prompt 3: query peptide ────────────────────────────────────────────────
  cat("\n[3/4] Enter your QUERY PEPTIDE SEQUENCE (single-letter amino acid codes):\n")
  query <- prompt_text_input(
    "Query peptide",
    validator  = function(x) nchar(x) > 0 && grepl("^[A-Za-z]+$", x),
    error_msg  = "Please enter a valid amino acid sequence (letters only, no spaces)."
  )
  query <- toupper(query)
  cat(sprintf('  -> Query: %s  (%d aa)\n', query, nchar(query)))

  # ── Prompt 4: top N ────────────────────────────────────────────────────────
  cat("\n[4/4] How many TOP RESULTS do you want per ranking? (default: 10)\n")
  top_n_input <- prompt_text_input(
    "Top N",
    validator = function(x) !is.na(suppressWarnings(as.integer(x))) && as.integer(x) > 0,
    error_msg = "Please enter a positive whole number."
  )
  top_n <- as.integer(top_n_input)
  cat(sprintf("  -> Returning top %d results per ranking.\n", top_n))

  # ── Filter dataset ─────────────────────────────────────────────────────────
  cat("\n  Filtering candidates (", norm_col, " > 0) ...\n", sep = "")

  df_filt <- df[
    !is.na(df[[norm_col]]) & df[[norm_col]] > 0 &
    !is.na(df[[pep_col]])  & nchar(as.character(df[[pep_col]])) > 0 &
    toupper(as.character(df[[pep_col]])) != query, ]

  if (nrow(df_filt) == 0) stop("No peptides pass the intensity > 0 filter.")
  cat(sprintf("  %d candidate peptides found.\n", nrow(df_filt)))

  # ── Score all candidates ───────────────────────────────────────────────────
  cat("  Scoring candidates ...\n")
  query_aas  <- strsplit(query, "")[[1]]
  candidates <- toupper(as.character(df_filt[[pep_col]]))

  all_scores <- lapply(candidates, function(cand) compute_scores(query_aas, cand))

  df_filt$norm_similarity <- sapply(all_scores, `[[`, "norm_similarity")
  df_filt$raw_similarity  <- sapply(all_scores, `[[`, "raw_similarity")
  df_filt$norm_identity   <- sapply(all_scores, `[[`, "norm_identity")
  df_filt$n_identical     <- sapply(all_scores, `[[`, "n_identical")
  df_filt$len_diff        <- sapply(all_scores, `[[`, "len_diff")
  df_filt$cand_length     <- sapply(all_scores, `[[`, "cand_length")
  df_filt$position_scores <- sapply(all_scores, `[[`, "position_scores")
  df_filt$query_length    <- nchar(query)

  # ── Build output columns: ONLY pep_col + norm_col + annotation + scores ────
  # Deliberately exclude all other intensity columns
  annot_cols <- c("Gene", "Protein.Description", "Protein.ID", "Entry.Name",
                  "Protein", "Mapped.Genes")
  annot_cols <- annot_cols[annot_cols %in% colnames(df_filt)]

  base_cols <- unique(c(
    pep_col,
    annot_cols,
    "query_length", "cand_length", "len_diff",
    "n_identical", "position_scores",
    norm_col          # only the user-specified normal column
  ))
  base_cols <- base_cols[base_cols %in% colnames(df_filt)]

  # ── Top SIMILARITY table ───────────────────────────────────────────────────
  sim_sorted <- df_filt[order(-df_filt$norm_similarity, -df_filt$raw_similarity), ]
  top_sim    <- head(sim_sorted[, c(base_cols, "norm_similarity",
                                     "raw_similarity"), drop = FALSE], top_n)
  rownames(top_sim) <- NULL

  # ── Top IDENTITY table ─────────────────────────────────────────────────────
  id_sorted <- df_filt[order(-df_filt$norm_identity, -df_filt$n_identical), ]
  top_id    <- head(id_sorted[, c(base_cols, "norm_identity",
                                   "n_identical"), drop = FALSE], top_n)
  rownames(top_id) <- NULL

  # ── Print results ──────────────────────────────────────────────────────────
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════╗\n")
  cat("║                      Results                            ║\n")
  cat("╚══════════════════════════════════════════════════════════╝\n")
  cat(sprintf("  Query          : %s  (%d aa)\n", query, nchar(query)))
  cat(sprintf("  Intensity col  : %s\n", norm_col))
  cat(sprintf("  Candidates     : %d\n\n", nrow(df_filt)))
  cat("  Scoring:\n")
  cat("    norm_similarity : identical(2pts) + conservative(1pt) - length penalty, /[0,1]\n")
  cat("    norm_identity   : exact matches only                  - length penalty, /[0,1]\n")
  cat("    position_scores : per-position score string (2|1|0), left to right\n")
  cat("\n  AA Groups for conservative substitution:\n")
  for (nm in names(AA_GROUPS)) {
    cat(sprintf("    %-20s: %s\n", nm, paste(AA_GROUPS[[nm]], collapse = " ")))
  }

  cat(sprintf("\n── TOP %d by SEQUENCE SIMILARITY ────────────────────────────\n", top_n))
  print(top_sim)

  cat(sprintf("\n── TOP %d by SEQUENCE IDENTITY ──────────────────────────────\n", top_n))
  print(top_id)

  invisible(list(top_similarity = top_sim,
                 top_identity   = top_id))
}

# ── 6. Alignment printer ──────────────────────────────────────────────────────
print_alignment <- function(query, candidate) {
  q       <- strsplit(toupper(query),     "")[[1]]
  cand    <- strsplit(toupper(candidate), "")[[1]]
  overlap <- min(length(q), length(cand))

  sym <- sapply(seq_len(overlap), function(i) {
    s <- score_position(q[i], cand[i])
    if (s == 2L) "|" else if (s == 1L) "+" else "."
  })

  n_id    <- sum(sym == "|")
  n_cons  <- sum(sym == "+")
  pct_id  <- round(100 * n_id / length(q), 1)
  pct_sim <- round(100 * (n_id + n_cons) / length(q), 1)

  cat("\nAlignment  (| identical  + conservative  . non-conservative)\n")
  cat("Query    :", paste(q[1:overlap],    collapse = " "), "\n")
  cat("         :", paste(sym,             collapse = " "), "\n")
  cat("Candidate:", paste(cand[1:overlap], collapse = " "), "\n")
  cat(sprintf("\n  Identity   : %d/%d (%.1f%%)\n", n_id, length(q), pct_id))
  cat(sprintf("  Similarity : %d/%d (%.1f%%)  [identical + conservative]\n\n",
              n_id + n_cons, length(q), pct_sim))
}

# ─────────────────────────────────────────────────────────────────────────────
# USAGE
# ─────────────────────────────────────────────────────────────────────────────
# source("peptide_similarity_search.R")
# results <- run_peptide_search(LMD)
#
# # You will be prompted to enter:
# #   [1] Peptide sequence column  e.g.  Peptide.Sequence
# #   [2] Normal intensity column  e.g.  Norm1.1
# #   [3] Query peptide            e.g.  AAAALVLKA
# #   [4] Top N                    e.g.  10
#
# # Access results:
# results$top_similarity
# results$top_identity
#
# # View alignment:
# print_alignment("AAAALVLKA", results$top_similarity$Peptide.Sequence[1])
# ─────────────────────────────────────────────────────────────────────────────
