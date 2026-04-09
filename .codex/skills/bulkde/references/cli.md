# bulkde CLI Reference

Use this file only when the task needs exact defaults, flags, or output filenames.

## Program shape

- Go CLI wrapper around an embedded R script
- Entry pattern: `bulkde run ...`
- Source tree in this environment: `/Volumes/DataCenter_01/msy0316/P25112107_2_bulk_result/02.quant/bulkde`

## Required flags

- `--counts`: count matrix path
- `--out`: output directory

## Important defaults

- `--group-from auto`
- `--gene-id-col ID`
- `--annot-prefix gene_`
- `--gene-name-col gene_name`
- `--meta-sample-col sample`
- `--group-col group`
- `--marker FOLH1`
- `--marker-field gene_name`
- `--marker-threshold 0`
- `--marker-op gt`
- `--min-count 10`
- `--min-samples 3`
- `--var-drop-quantile 0.25`
- `--methods limma,deseq2,edger`
- `--fdr 0.05`
- `--lfc 1`

## Grouping modes

### Marker mode

Use when there is no metadata table or when the user explicitly wants marker-high versus marker-low grouping.

- `--group-from marker`
- `--marker`
- `--marker-field gene_name|gene_id`
- `--marker-threshold`
- `--marker-op gt|ge|eq|ne|lt|le`

Default interpretation is `marker counts > threshold => POS`, otherwise `NEG`.

### Metadata mode

Use when the user has a sample annotation file.

- `--group-from meta`
- `--meta`
- `--meta-sample-col`
- `--group-col`
- `--case`
- `--control`
- `--covariates`

If the group column has exactly two levels, the tool can infer case/control when the user does not provide them. If there are more than two levels, the comparison must be narrowed explicitly.

## Common outputs

- `sample_meta.tsv`
- `sample_group.tsv`
- `run_config.tsv`
- `gene_count.filtered.tsv`
- `filter_summary.tsv`
- `de_summary.tsv`
- `session_info.txt`
- `limma_all.tsv`, `limma_sig.tsv`
- `deseq2_all.tsv`, `deseq2_sig.tsv`
- `edger_all.tsv`, `edger_sig.tsv`

## Runtime dependencies

- `R`
- `rvx` from `rs-reborn`, or legacy `rs`

`bulkde` prefers `rvx` because some systems expose a different unrelated `/usr/bin/rs`.

## Common failure patterns

- `--counts/--count is required`: required input path missing.
- `--out is required`: output directory missing.
- `No sample columns detected`: annotation columns were mistaken for sample columns, or count columns are not numeric.
- `Marker row not found`: wrong marker value or wrong `--marker-field`.
- `Meta missing samples`: sample names do not match between counts and metadata.
- `Expression filter removed all genes`: thresholds are too strict for the matrix.
