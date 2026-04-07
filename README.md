# bulkde

A small Go CLI wrapper that runs bulk RNA-seq differential expression via an embedded R script (limma-voom, DESeq2, edgeR).

It performs:

- Grouping:
  - marker mode: `marker counts > threshold` (or chosen op) => `POS`, else `NEG`
  - meta mode: uses `meta` table `group` column and compares `case - control` (supports covariates)
- Prefiltering: keep genes with `counts >= min-count` in at least `min-samples` samples
- Low-variance removal: drop the lowest `var-drop-quantile` by `var(log2(CPM+1))`
- DE methods (same filtered genes): `limma-voom`, `DESeq2`, `edgeR (QL)`

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

- Requires `R` and the `rs` CLI from `rs-reborn` on PATH. `bulkde` executes `rs run ...` to resolve R packages into an isolated cache, then runs the embedded R script.
- `--cache-dir` controls the rs-reborn cache root (default: `<counts_dir>/r_libs`). `--r-lib` is a backwards-compatible alias. You can override the `rs` binary path via `--rs-path`.
- Some systems already ship a different `/usr/bin/rs` tool (not rs-reborn). `bulkde` will validate `rs run --help` is available, otherwise it will error and ask you to pass `--rs-path`.
- Backwards compatible alias: `--count` works as alias of `--counts`.
