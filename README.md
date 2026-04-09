# bulkde

A small Go CLI wrapper that runs bulk RNA-seq differential expression via an embedded R script (limma-voom, DESeq2, edgeR), with an optional `--chip-mode` for preprocessed microarray/chip expression matrices.

It performs:

- Grouping:
  - marker mode: `marker counts > threshold` (or chosen op) => `POS`, else `NEG`
  - meta mode: uses `meta` table `group` column and compares `case - control` (supports covariates)
- Prefiltering: keep genes with `counts >= min-count` in at least `min-samples` samples
- Low-variance removal: drop the lowest `var-drop-quantile` by `var(log2(CPM+1))`
- DE methods (same filtered genes): `limma-voom`, `DESeq2`, `edgeR (QL)`
- Chip mode: `--chip-mode` treats the input as a preprocessed chip expression matrix, forces `limma` only, and skips RNA-seq-specific count filtering

## Build

```bash
cd /Volumes/DataCenter_01/msy0316/P25112107_2_bulk_result/02.quant/bulkde
GOCACHE=/tmp/go-build go build -o bulkde
```

## Run

Example (FOLH1 POS vs NEG):

```bash
./bulkde run \
  --counts /Volumes/DataCenter_01/msy0316/P25112107_2_bulk_result/02.quant/gene_count.txt \
  --out    /path/to/out_dir \
  --group-from marker \
  --marker FOLH1 \
  --no-install
```

Chip / microarray example:

```bash
./bulkde run \
  --counts /path/to/chip_expression.tsv \
  --out /path/to/out_dir \
  --group-from meta \
  --meta /path/to/meta.tsv \
  --case Tumor \
  --control Normal \
  --chip-mode \
  --no-install
```

Outputs under `--out`:

- `sample_meta.tsv` (+ compat `sample_group.tsv`)
- `run_config.tsv`
- `gene_count.filtered.tsv`
- `filter_summary.tsv`
- `limma_all.tsv`, `limma_sig.tsv`
- `deseq2_all.tsv`, `deseq2_sig.tsv`
- `edger_all.tsv`, `edger_sig.tsv`
- `de_summary.tsv`
- `session_info.txt`

## Notes

- Requires `R` and the runner CLI from `rs-reborn` on PATH. The upstream binary was renamed from `rs` to `rvx`; `bulkde` will prefer `rvx` and can fall back to `rs` if present.
- `bulkde` executes `rvx run ...` to resolve R packages into an isolated cache, then runs the embedded R script.
- `--cache-dir` controls the rs-reborn cache root (default: `<counts_dir>/r_libs`). `--r-lib` is a backwards-compatible alias.
- You can override the runner binary path via `--rvx-path` (`--rs-path` is a compat alias).
- Some systems ship a different `/usr/bin/rs` tool (not rs-reborn). Prefer using `rvx` to avoid name conflicts.
- Backwards compatible alias: `--count` works as alias of `--counts`.
- `--chip-mode` is intended for already preprocessed chip/microarray expression matrices and runs `limma` only.
- `bulkde` does not perform chip probe annotation, probe-to-gene collapsing, or duplicated-gene deduplication for chip inputs. Please complete probe annotation and gene deduplication before running `--chip-mode`.
