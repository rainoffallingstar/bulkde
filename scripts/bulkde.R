#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: bulkde.R <counts.tsv> <out_dir>\n", file = stderr())
  quit(status = 2)
}

counts_file <- args[[1]]
out_dir <- args[[2]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

getenv1 <- function(key, default = "") {
  v <- Sys.getenv(key, unset = default)
  if (is.na(v) || !nzchar(v)) default else v
}

as_bool <- function(x) tolower(x) %in% c("1", "true", "t", "yes", "y")

# rs-reborn manages the library location and injects it via .libPaths().
local_lib <- .libPaths()[[1]]
no_install <- as_bool(getenv1("BULKDE_NO_INSTALL", "0"))

# Config (from env)
gene_id_col <- getenv1("BULKDE_GENE_ID_COL", "ID")
annot_prefix <- getenv1("BULKDE_ANNOT_PREFIX", "gene_")
annot_cols_raw <- getenv1("BULKDE_ANNOT_COLS", "")
gene_name_col <- getenv1("BULKDE_GENE_NAME_COL", "gene_name")

group_from <- tolower(getenv1("BULKDE_GROUP_FROM", "auto")) # auto|meta|marker
meta_file <- getenv1("BULKDE_META", "")
meta_sample_col <- getenv1("BULKDE_META_SAMPLE_COL", "sample")
group_col <- getenv1("BULKDE_GROUP_COL", "group")
case_label <- getenv1("BULKDE_CASE", "")
control_label <- getenv1("BULKDE_CONTROL", "")
covariates_raw <- getenv1("BULKDE_COVARIATES", "")

marker_gene <- getenv1("BULKDE_MARKER", "FOLH1")
marker_field <- tolower(getenv1("BULKDE_MARKER_FIELD", "gene_name")) # gene_name|gene_id
marker_threshold <- as.numeric(getenv1("BULKDE_MARKER_THRESHOLD", "0"))
marker_op <- tolower(getenv1("BULKDE_MARKER_OP", "gt")) # gt|ge|eq|ne|lt|le

min_count <- as.integer(getenv1("BULKDE_MIN_COUNT", "10"))
min_samples <- as.integer(getenv1("BULKDE_MIN_SAMPLES", "3"))
var_drop_quantile <- as.numeric(getenv1("BULKDE_VAR_DROP_QUANTILE", "0.25"))

methods_raw <- tolower(getenv1("BULKDE_METHODS", "limma,deseq2,edger"))
fdr_threshold <- as.numeric(getenv1("BULKDE_FDR", "0.05"))
lfc_threshold <- as.numeric(getenv1("BULKDE_LFC", "1"))

split_csv <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(character())
  parts <- unlist(strsplit(x, ",", fixed = TRUE))
  parts <- trimws(parts)
  parts[parts != "" & !is.na(parts)]
}

methods <- unique(split_csv(methods_raw))
methods <- methods[methods %in% c("limma", "deseq2", "edger")]
if (length(methods) == 0) stop("No valid methods selected. Allowed: limma,deseq2,edger")

suppressPackageStartupMessages({
  library(edgeR)
  if ("limma" %in% methods) library(limma)
  if ("deseq2" %in% methods) library(DESeq2)
})

covariates <- split_csv(covariates_raw)

annot_cols_user <- split_csv(annot_cols_raw)
# Sensible defaults for common count-matrix formats (e.g. featureCounts output).
annot_cols_default <- c(
  # featureCounts
  "Chr", "Start", "End", "Strand", "Length",
  # StringTie (some matrices include extra non-sample columns)
  "gene_name", "gene", "gene_id", "gene_short_name", "gene_shortname", "symbol", "description",
  "ref_gene_id", "ref_gene_name", "reference_id", "reference", "transcript_id",
  "locus", "chrom", "chromosome",
  # RSEM-like columns that should never be treated as sample columns
  "effective_length", "expected_count", "tpm", "fpkm", "isopct"
)
annot_cols_explicit <- unique(c(annot_cols_default, annot_cols_user))

detect_sep <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, n = 50, warn = FALSE)
  line <- ""
  for (l in lines) {
    if (!startsWith(l, "#") && nzchar(trimws(l))) {
      line <- l
      break
    }
  }
  if (!nzchar(line)) return("\t")
  if (grepl("\t", line, fixed = TRUE)) return("\t")
  if (grepl(",", line, fixed = TRUE)) return(",")
  "\t"
}

read_meta <- function(path) {
  if (!nzchar(path)) return(NULL)
  if (!file.exists(path)) stop("Meta file not found: ", path)
  first <- readLines(path, n = 1, warn = FALSE)
  sep <- if (length(first) == 1 && grepl(",", first, fixed = TRUE) && !grepl("\t", first, fixed = TRUE)) "," else "\t"
  read.table(path, header = TRUE, sep = sep, quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
}

sep_counts <- detect_sep(counts_file)
con_counts <- if (grepl("\\.gz$", counts_file, ignore.case = TRUE)) gzfile(counts_file, open = "rt") else file(counts_file, open = "rt")
counts_df <- read.table(con_counts, header = TRUE, sep = sep_counts, quote = "", comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)
close(con_counts)
if (nrow(counts_df) == 0 || ncol(counts_df) < 2) stop("Counts file looks empty or malformed: ", counts_file)

if (!(gene_id_col %in% colnames(counts_df))) {
  gene_id_col <- colnames(counts_df)[[1]]
}

is_annot <- function(col) {
  if (!nzchar(annot_prefix)) return(FALSE)
  startsWith(col, annot_prefix)
}

all_cols <- colnames(counts_df)
# Case-insensitive match for explicit annotation column names
all_cols_lower <- tolower(all_cols)
explicit_lower <- tolower(annot_cols_explicit)
explicit_hits <- all_cols[all_cols_lower %in% explicit_lower]

annot_cols <- unique(c(all_cols[is_annot(all_cols)], explicit_hits))
candidate_cols <- setdiff(all_cols, c(gene_id_col, annot_cols))

is_numeric_col <- function(v) {
  if (is.numeric(v) || is.integer(v)) return(TRUE)
  x <- as.character(v)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(FALSE)
  suppressWarnings(num <- as.numeric(x))
  mean(!is.na(num)) >= 0.9
}

sample_cols <- candidate_cols[vapply(candidate_cols, function(cn) is_numeric_col(counts_df[[cn]]), logical(1))]
if (length(sample_cols) == 0) stop("No sample columns detected in counts after excluding gene_id and annot columns.")

gene_ids <- as.character(counts_df[[gene_id_col]])
if (any(is.na(gene_ids) | gene_ids == "")) stop("Empty gene_id values detected in column: ", gene_id_col)

# HTSeq-count summary rows start with "__" (not real genes); drop them if present.
is_htseq_summary <- startsWith(gene_ids, "__")
if (any(is_htseq_summary, na.rm = TRUE)) {
  counts_df <- counts_df[!is_htseq_summary, , drop = FALSE]
  gene_ids <- gene_ids[!is_htseq_summary]
}

counts_mat <- as.matrix(counts_df[, sample_cols, drop = FALSE])
storage.mode(counts_mat) <- "integer"
rownames(counts_mat) <- gene_ids

gene_annot <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE, check.names = FALSE)
if (length(annot_cols) > 0) {
  gene_annot <- cbind(gene_annot, counts_df[, annot_cols, drop = FALSE])
}

get_marker_row_idx <- function() {
  if (marker_field == "gene_id") {
    ix <- which(gene_ids == marker_gene)
    if (length(ix) != 1) return(integer())
    return(ix)
  }
  if (marker_field == "gene_name") {
    if (!(gene_name_col %in% colnames(counts_df))) return(integer())
    ix <- which(as.character(counts_df[[gene_name_col]]) == marker_gene)
    if (length(ix) != 1) return(integer())
    return(ix)
  }
  integer()
}

marker_idx <- get_marker_row_idx()
marker_counts <- rep(NA_real_, length(sample_cols))
if (length(marker_idx) == 1) {
  marker_counts <- as.numeric(counts_df[marker_idx, sample_cols, drop = TRUE])
}

apply_marker_op <- function(x) {
  if (marker_op == "gt") return(x > marker_threshold)
  if (marker_op == "ge") return(x >= marker_threshold)
  if (marker_op == "eq") return(x == marker_threshold)
  if (marker_op == "ne") return(x != marker_threshold)
  if (marker_op == "lt") return(x < marker_threshold)
  if (marker_op == "le") return(x <= marker_threshold)
  stop("Unsupported marker-op: ", marker_op)
}

if (group_from == "auto") {
  group_from <- if (nzchar(meta_file)) "meta" else "marker"
}
if (!group_from %in% c("meta", "marker")) stop("Invalid group-from: ", group_from)

meta_df <- NULL
group_used <- NULL
case_used <- NULL
control_used <- NULL

if (group_from == "marker") {
  if (length(marker_idx) != 1) {
    stop("Marker row not found (or not unique) for marker=", marker_gene, " marker-field=", marker_field, ". Needed for marker grouping.")
  }
  pos <- apply_marker_op(marker_counts)
  group_used <- ifelse(pos, "POS", "NEG")
  case_used <- "POS"
  control_used <- "NEG"
  meta_df <- data.frame(
    sample = sample_cols,
    group = group_used,
    stringsAsFactors = FALSE
  )
} else {
  meta_in <- read_meta(meta_file)
  if (is.null(meta_in) || nrow(meta_in) == 0) stop("Meta file is empty: ", meta_file)
  if (!(meta_sample_col %in% colnames(meta_in))) stop("Missing meta sample column: ", meta_sample_col)
  if (!(group_col %in% colnames(meta_in))) stop("Missing meta group column: ", group_col)

  meta_in[[meta_sample_col]] <- as.character(meta_in[[meta_sample_col]])
  meta_in <- meta_in[!is.na(meta_in[[meta_sample_col]]) & meta_in[[meta_sample_col]] != "", , drop = FALSE]
  meta_in <- meta_in[!duplicated(meta_in[[meta_sample_col]]), , drop = FALSE]

  missing <- setdiff(sample_cols, meta_in[[meta_sample_col]])
  if (length(missing) > 0) stop("Meta missing samples: ", paste(missing, collapse = ", "))

  meta_df <- meta_in[match(sample_cols, meta_in[[meta_sample_col]]), , drop = FALSE]
  colnames(meta_df)[colnames(meta_df) == meta_sample_col] <- "sample"
  meta_df$group <- as.character(meta_df[[group_col]])

  # Determine case/control if not provided and there are exactly 2 groups.
  grp_levels <- unique(meta_df$group)
  if (!nzchar(control_label) || !nzchar(case_label)) {
    if (length(grp_levels) != 2) {
      stop("Meta grouping requires --case and --control when group has != 2 levels. Found: ", paste(grp_levels, collapse = ", "))
    }
    control_used <- grp_levels[[1]]
    case_used <- grp_levels[[2]]
  } else {
    control_used <- control_label
    case_used <- case_label
  }

  keep <- meta_df$group %in% c(control_used, case_used)
  if (!all(keep)) {
    dropped <- unique(meta_df$group[!keep])
    meta_df <- meta_df[keep, , drop = FALSE]
    counts_mat <- counts_mat[, meta_df$sample, drop = FALSE]
    sample_cols <- meta_df$sample
    marker_counts <- marker_counts[match(sample_cols, meta_df$sample)]
    message("Dropping samples not in case/control groups: ", paste(dropped, collapse = ", "))
  }

  if (!all(c(control_used, case_used) %in% unique(meta_df$group))) {
    stop("Meta group does not contain both case and control. case=", case_used, " control=", control_used)
  }
}

# Attach marker_count if available (even for meta grouping)
meta_df$marker_gene <- marker_gene
meta_df$marker_count <- marker_counts[match(meta_df$sample, sample_cols)]

# Attach covariates if requested
if (length(covariates) > 0) {
  for (cv in covariates) {
    if (!(cv %in% colnames(meta_df))) stop("Missing covariate in meta: ", cv)
  }
}

# Build analysis frame with only needed columns, and coerce types.
analysis_df <- meta_df
analysis_df$group <- factor(analysis_df$group, levels = c(control_used, case_used))
for (cv in covariates) {
  v <- analysis_df[[cv]]
  if (is.character(v)) {
    # Detect numeric
    suppressWarnings(num <- as.numeric(v))
    if (!all(is.na(num) & !is.na(v)) && sum(!is.na(num)) == sum(!is.na(v))) {
      analysis_df[[cv]] <- num
    } else {
      analysis_df[[cv]] <- factor(v)
    }
  }
}

pos_n <- sum(analysis_df$group == case_used)
neg_n <- sum(analysis_df$group == control_used)
if (pos_n < 1 || neg_n < 1) stop("Need at least 1 sample in each group. case=", pos_n, " control=", neg_n)

# Outputs: sample_meta + backwards-compatible sample_group
sample_meta_out <- data.frame(
  sample = analysis_df$sample,
  group = as.character(analysis_df$group),
  role = ifelse(as.character(analysis_df$group) == case_used, "case", "control"),
  case = case_used,
  control = control_used,
  marker_gene = analysis_df$marker_gene,
  marker_count = analysis_df$marker_count,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (cv in covariates) sample_meta_out[[cv]] <- analysis_df[[cv]]
write.table(sample_meta_out, file = file.path(out_dir, "sample_meta.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

sample_group_compat <- sample_meta_out[, c("sample", "marker_count", "group")]
write.table(sample_group_compat, file = file.path(out_dir, "sample_group.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Record run_config
cfg <- data.frame(
  key = c(
    "counts_file", "out_dir", "gene_id_col", "annot_prefix", "gene_name_col",
    "group_from", "meta_file", "meta_sample_col", "group_col",
    "case", "control", "covariates",
    "marker_gene", "marker_field", "marker_op", "marker_threshold",
    "min_count", "min_samples", "var_drop_quantile",
    "methods", "fdr_threshold", "lfc_threshold",
    "r_lib", "no_install"
  ),
  value = c(
    normalizePath(counts_file, winslash = "/", mustWork = FALSE),
    normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    gene_id_col, annot_prefix, gene_name_col,
    group_from, if (nzchar(meta_file)) normalizePath(meta_file, winslash = "/", mustWork = FALSE) else "",
    meta_sample_col, group_col,
    case_used, control_used, paste(covariates, collapse = ","),
    marker_gene, marker_field, marker_op, as.character(marker_threshold),
    as.character(min_count), as.character(min_samples), as.character(var_drop_quantile),
    paste(methods, collapse = ","), as.character(fdr_threshold), as.character(lfc_threshold),
    normalizePath(local_lib, winslash = "/", mustWork = FALSE),
    ifelse(no_install, "1", "0")
  ),
  stringsAsFactors = FALSE
)
write.table(cfg, file = file.path(out_dir, "run_config.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Filtering
keep_expression <- rowSums(counts_mat >= min_count) >= min_samples
counts_expr <- counts_mat[keep_expression, , drop = FALSE]
annot_expr <- gene_annot[keep_expression, , drop = FALSE]
if (nrow(counts_expr) == 0) stop("Expression filter removed all genes; please revisit thresholds.")

counts_filtered <- counts_expr
annot_filtered <- annot_expr
if (var_drop_quantile > 0) {
  log_cpm <- cpm(counts_expr, log = TRUE, prior.count = 1)
  gene_var <- apply(log_cpm, 1, var)
  var_cutoff <- unname(stats::quantile(gene_var, probs = var_drop_quantile, na.rm = TRUE, type = 7))
  keep_variance <- gene_var > var_cutoff
  if (!any(keep_variance)) stop("Variance filter removed all genes; please revisit thresholds.")
  counts_filtered <- counts_expr[keep_variance, , drop = FALSE]
  annot_filtered <- annot_expr[keep_variance, , drop = FALSE]
}

filtered_export <- data.frame(gene_id = rownames(counts_filtered), counts_filtered, stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(annot_filtered) > 1) {
  # annot_filtered includes gene_id + annot cols
  filtered_export <- cbind(filtered_export, annot_filtered[, setdiff(colnames(annot_filtered), "gene_id"), drop = FALSE])
}
write.table(filtered_export, file = file.path(out_dir, "gene_count.filtered.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

filter_summary <- data.frame(
  metric = c("raw_gene_count", "after_expression_filter", "after_variance_filter"),
  value = c(nrow(counts_df), nrow(counts_expr), nrow(counts_filtered)),
  stringsAsFactors = FALSE
)
write.table(filter_summary, file = file.path(out_dir, "filter_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

merge_with_annot <- function(res_df) {
  merge(annot_filtered, res_df, by = "gene_id", all.y = TRUE, sort = FALSE)
}

add_comparison_cols <- function(df) {
  df$case <- case_used
  df$control <- control_used
  df$comparison <- paste0(case_used, "_vs_", control_used)
  df
}

write_sig_table <- function(df, file_name, padj_col = "FDR") {
  sig <- df[!is.na(df[[padj_col]]) & df[[padj_col]] < fdr_threshold & abs(df$log2FC) >= lfc_threshold, , drop = FALSE]
  write.table(sig, file = file.path(out_dir, file_name), sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(sig)
}

summarise_sig <- function(df, method, padj_col) {
  sig <- df[!is.na(df[[padj_col]]) & df[[padj_col]] < fdr_threshold & abs(df$log2FC) >= lfc_threshold, , drop = FALSE]
  data.frame(
    method = method,
    all_genes = nrow(df),
    sig_genes = nrow(sig),
    up = sum(sig$log2FC > 0, na.rm = TRUE),
    down = sum(sig$log2FC < 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# Design matrix/formula
design_terms <- c(covariates, "group")
design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))

design <- model.matrix(design_formula, data = analysis_df)
rownames(design) <- analysis_df$sample
case_coef_name <- paste0("group", make.names(case_used))
coef_idx <- which(colnames(design) == case_coef_name)
if (length(coef_idx) != 1) {
  stop("Could not find group coefficient in design matrix. Expected column: ", case_coef_name, " got: ", paste(colnames(design), collapse = ", "))
}
contrast <- rep(0, ncol(design))
contrast[coef_idx] <- 1

# Persist design and contrast for reproducibility and debugging; also print key info to logs.
design_out <- data.frame(sample = analysis_df$sample, design, stringsAsFactors = FALSE, check.names = FALSE)
write.table(design_out, file = file.path(out_dir, "design_matrix.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
contrast_out <- data.frame(term = colnames(design), weight = as.numeric(contrast), stringsAsFactors = FALSE, check.names = FALSE)
write.table(contrast_out, file = file.path(out_dir, "contrast.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

message("Comparison (log2FC direction): ", case_used, " - ", control_used)
message("Interpretation: log2FC > 0 means higher expression in case (", case_used, ") compared to control (", control_used, ").")
message("Group factor levels (baseline first): ", paste(levels(analysis_df$group), collapse = " < "))
message("Design formula: ", paste(deparse(design_formula), collapse = " "))
message("Design matrix dim: ", nrow(design), " x ", ncol(design))
message("Design columns: ", paste(colnames(design), collapse = ", "))
message("Contrast column: ", colnames(design)[coef_idx])
message("Contrast weights: ", paste(sprintf("%s:%g", colnames(design), as.numeric(contrast)), collapse = ", "))
message("Design matrix preview (first 6 rows):\n", paste(capture.output(utils::head(design)), collapse = "\n"))

deseq_out <- NULL
edger_out <- NULL
limma_out <- NULL

dge <- DGEList(counts = counts_filtered)
dge <- calcNormFactors(dge)

if ("limma" %in% methods) {
  voom_obj <- voom(dge, design = design, plot = FALSE)
  fit <- lmFit(voom_obj, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit)
  limma_tbl <- topTable(fit, coef = 1, number = Inf, sort.by = "P")
  limma_tbl <- data.frame(gene_id = rownames(limma_tbl), limma_tbl, row.names = NULL, check.names = FALSE)
  limma_tbl <- transform(limma_tbl, log2FC = logFC, pvalue = P.Value, padj = adj.P.Val)
  limma_out <- limma_tbl[, c("gene_id", "log2FC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "pvalue", "padj")]
  colnames(limma_out)[colnames(limma_out) == "adj.P.Val"] <- "FDR"
  limma_out <- merge_with_annot(limma_out)
  limma_out <- add_comparison_cols(limma_out)
  write.table(limma_out, file = file.path(out_dir, "limma_all.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write_sig_table(limma_out, "limma_sig.tsv", padj_col = "FDR")
}

if ("deseq2" %in% methods) {
  # Build colData with proper rownames
  coldata <- analysis_df[, c(covariates, "group"), drop = FALSE]
  rownames(coldata) <- analysis_df$sample
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered[, rownames(coldata), drop = FALSE],
    colData = coldata,
    design = design_formula
  )
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  deseq_res <- results(dds, contrast = c("group", case_used, control_used))
  deseq_tbl <- as.data.frame(deseq_res)
  deseq_tbl <- data.frame(gene_id = rownames(deseq_tbl), deseq_tbl, row.names = NULL, check.names = FALSE)
  deseq_tbl <- transform(deseq_tbl, log2FC = log2FoldChange, pvalue = pvalue, padj = padj)
  deseq_out <- deseq_tbl[, c("gene_id", "baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj")]
  colnames(deseq_out)[colnames(deseq_out) == "padj"] <- "FDR"
  deseq_out <- merge_with_annot(deseq_out)
  deseq_out <- add_comparison_cols(deseq_out)
  write.table(deseq_out, file = file.path(out_dir, "deseq2_all.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write_sig_table(deseq_out, "deseq2_sig.tsv", padj_col = "FDR")
}

if ("edger" %in% methods) {
  dge_edger <- dge
  dge_edger <- estimateDisp(dge_edger, design)
  qlf_fit <- glmQLFit(dge_edger, design)
  qlf <- glmQLFTest(qlf_fit, contrast = contrast)
  edger_tbl <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  edger_tbl <- data.frame(gene_id = rownames(edger_tbl), edger_tbl, row.names = NULL, check.names = FALSE)
  edger_tbl <- transform(edger_tbl, log2FC = logFC, pvalue = PValue, padj = FDR)
  edger_out <- edger_tbl[, c("gene_id", "log2FC", "logCPM", "F", "PValue", "FDR", "pvalue", "padj")]
  edger_out <- merge_with_annot(edger_out)
  edger_out <- add_comparison_cols(edger_out)
  write.table(edger_out, file = file.path(out_dir, "edger_all.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  write_sig_table(edger_out, "edger_sig.tsv", padj_col = "FDR")
}

# Summary across executed methods
sum_rows <- list()
if (!is.null(limma_out)) sum_rows[["limma"]] <- summarise_sig(limma_out, "limma", "FDR")
if (!is.null(deseq_out)) sum_rows[["DESeq2"]] <- summarise_sig(deseq_out, "DESeq2", "FDR")
if (!is.null(edger_out)) sum_rows[["edgeR"]] <- summarise_sig(edger_out, "edgeR", "FDR")
de_summary <- do.call(rbind, sum_rows)
write.table(de_summary, file = file.path(out_dir, "de_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

session_file <- file.path(out_dir, "session_info.txt")
sink(session_file)
cat("Analysis date:", format(Sys.time(), tz = Sys.timezone()), "\n")
cat("Counts file:", normalizePath(counts_file, winslash = "/", mustWork = FALSE), "\n")
cat("Output dir:", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n")
cat("Grouping: group_from=", group_from, " case=", case_used, " control=", control_used, "\n", sep = "")
cat("Design:", deparse(design_formula), "\n")
cat("Design matrix file:", normalizePath(file.path(out_dir, "design_matrix.tsv"), winslash = "/", mustWork = FALSE), "\n")
cat("Contrast file:", normalizePath(file.path(out_dir, "contrast.tsv"), winslash = "/", mustWork = FALSE), "\n")
cat("Samples: n=", nrow(analysis_df), " case=", pos_n, " control=", neg_n, "\n", sep = "")
cat("Marker:", marker_gene, " field=", marker_field, " op=", marker_op, " threshold=", marker_threshold, " row_found=", length(marker_idx) == 1, "\n", sep = "")
cat("Filtering: min_count=", min_count, " min_samples=", min_samples, " var_drop_quantile=", var_drop_quantile, "\n", sep = "")
cat("Methods:", paste(methods, collapse = ","), "\n")
cat("Sig: FDR<", fdr_threshold, " |log2FC|>=", lfc_threshold, "\n\n", sep = "")
print(sessionInfo())
sink()

# Ensure required files exist
required_files <- c(
  "sample_meta.tsv", "sample_group.tsv", "run_config.tsv",
  "design_matrix.tsv", "contrast.tsv",
  "gene_count.filtered.tsv", "filter_summary.tsv", "de_summary.tsv", "session_info.txt"
)
for (m in methods) {
  if (m == "limma") required_files <- c(required_files, "limma_all.tsv", "limma_sig.tsv")
  if (m == "deseq2") required_files <- c(required_files, "deseq2_all.tsv", "deseq2_sig.tsv")
  if (m == "edger") required_files <- c(required_files, "edger_all.tsv", "edger_sig.tsv")
}
stopifnot(all(file.exists(file.path(out_dir, required_files))))

message("bulkde completed successfully.")
