---
name: "bulkde"
description: "Use when the user wants to run or interpret the local bulkde program for bulk RNA-seq differential expression, especially from a counts matrix plus optional sample metadata, marker-based POS/NEG grouping, or when they need limma, DESeq2, and edgeR outputs produced and summarized."
---

# Bulkde

Use this skill when the task is specifically about the local `bulkde` CLI in this workspace or a nearby checkout. This skill is for running bulk RNA-seq differential expression from a gene count matrix, choosing marker-based or metadata-based grouping, inspecting `bulkde` outputs, and troubleshooting common execution failures.

Do not use this skill for single-cell workflows, enrichment-only analysis, or generic R differential-expression code that does not rely on the `bulkde` program.

## Quick Start

1. Find the program first.
   Prefer an existing binary in the current workspace.
   Check `./bulkde`, `./bulkde/bulkde`, then the known source tree at `/Volumes/DataCenter_01/msy0316/P25112107_2_bulk_result/02.quant/bulkde`.
2. If no binary exists but the source tree is present, build it from the `bulkde/` directory with `GOCACHE=/tmp/go-build go build -o bulkde`.
3. Confirm runtime dependencies before running analysis.
   `bulkde` needs `R` plus the `rs-reborn` runner on `PATH`, usually `rvx` and sometimes legacy `rs`.
4. Choose a grouping mode.
   Use marker grouping when only a count matrix is available.
   Use metadata grouping when a sample table defines case/control and optional covariates.
5. Run `bulkde run ...` with an explicit output directory.
6. Read `de_summary.tsv` first, then inspect `*_sig.tsv`, `sample_meta.tsv`, and `filter_summary.tsv`.

## Workflow

### 1. Verify the inputs

Expect a count matrix with:

- one gene identifier column, default name `ID`
- optional annotation columns such as `gene_name`, `Chr`, `Start`, `End`, `Length`
- sample columns that are numeric counts

If metadata grouping is used, expect a TSV or CSV with:

- a sample column, default `sample`
- a group column, default `group`
- optional covariate columns such as `batch` or `sex`

When the user does not specify `gene-id-col`, `gene-name-col`, `group-col`, or `meta-sample-col`, rely on the program defaults unless the file headers clearly differ.

### 2. Choose grouping carefully

Use `--group-from marker` when the user wants POS vs NEG splitting from a marker gene. The default marker is `FOLH1`, default field is `gene_name`, and default split is `counts > 0 => POS`.

Use `--group-from meta` when the user has sample metadata and wants an explicit `case - control` comparison. If the metadata contains more than two groups, require `--case` and `--control`.

### 3. Prefer these run patterns

Marker-based example:

```bash
/absolute/path/to/bulkde run \
  --counts /absolute/path/to/gene_count.txt \
  --out /absolute/path/to/output_dir \
  --group-from marker \
  --marker FOLH1 \
  --marker-field gene_name \
  --no-install
```

Metadata-based example:

```bash
/absolute/path/to/bulkde run \
  --counts /absolute/path/to/gene_count.txt \
  --out /absolute/path/to/output_dir \
  --group-from meta \
  --meta /absolute/path/to/meta.tsv \
  --meta-sample-col sample \
  --group-col group \
  --case Tumor \
  --control Normal \
  --covariates batch,sex \
  --no-install
```

Use `--methods limma,deseq2,edger` unless the user wants a faster subset. Keep `--no-install` when the environment is already provisioned or package installation is likely to be restricted.

### 4. Review outputs in this order

1. `de_summary.tsv` for the top-level hit counts by method.
2. `sample_meta.tsv` to confirm group assignment and marker counts.
3. `filter_summary.tsv` to see whether filtering removed too many genes.
4. `limma_sig.tsv`, `deseq2_sig.tsv`, `edger_sig.tsv` for significant genes.
5. `*_all.tsv` when the user wants ranking beyond significance cutoffs.
6. `session_info.txt` for reproducibility and package debugging.

### 5. Troubleshoot with the shortest path

- If `rvx` or `rs` is missing, stop and surface that runtime dependency clearly.
- If no sample columns are detected, inspect the header and use `--annot-cols` or `--annot-prefix` so annotation columns are excluded from sample detection.
- If marker grouping fails, verify whether the marker should be matched through `gene_name` or `gene_id`.
- If metadata is missing samples, compare the count-matrix sample columns against the metadata sample column exactly.
- If all genes are filtered out, lower `--min-count`, lower `--min-samples`, or reduce `--var-drop-quantile`.
- If the user only needs one method, narrow `--methods` instead of debugging all three at once.

## Command Rules

- Always use absolute paths in the final execution command.
- Keep output directories separate per run so results are easy to inspect.
- Do not overwrite a user output directory blindly when it already contains prior results; prefer a sibling directory with a descriptive suffix.
- Summarize the comparison direction explicitly as `case - control` when reporting DE results.
- When interpreting results, mention which method table the claim comes from.

## Reference Map

Read [references/cli.md](references/cli.md) only when you need the full flag map, output inventory, or a quick reminder of defaults and common failure cases.
