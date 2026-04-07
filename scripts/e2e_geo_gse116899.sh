#!/usr/bin/env bash
set -euo pipefail

# End-to-end test using a real GEO series with processed featureCounts output.
# Dataset: GSE116899 (human neutrophils, featureCounts matrix with Geneid/Chr/Start/End/Strand/Length + samples).

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

GSE_ID="${BULKDE_E2E_GSE_ID:-GSE116899}"
URL="${BULKDE_E2E_URL:-https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116899/suppl/GSE116899_ra-neut-counts-EnsembIDs-GRCh37.p10.txt.gz}"
WORK_DIR="${BULKDE_E2E_WORK_DIR:-/tmp/bulkde-e2e-geo/${GSE_ID}}"
OUT_DIR="${BULKDE_E2E_OUT_DIR:-${WORK_DIR}/out}"
CACHE_DIR="${BULKDE_E2E_CACHE_DIR:-${WORK_DIR}/r_libs}"
N_GENES="${BULKDE_E2E_N_GENES:-10000}"   # subset for speed
METHODS="${BULKDE_E2E_METHODS:-limma}"   # keep fast by default
NO_INSTALL="${BULKDE_E2E_NO_INSTALL:-0}" # default to self-contained (may install R packages)

mkdir -p "${WORK_DIR}" "${OUT_DIR}"

RUNNER_BIN="${BULKDE_E2E_RUNNER_BIN:-}"
if [[ -z "${RUNNER_BIN}" ]]; then
  if command -v rvx >/dev/null 2>&1; then
    RUNNER_BIN="rvx"
  elif command -v rs >/dev/null 2>&1; then
    RUNNER_BIN="rs"
  fi
fi

if [[ -z "${RUNNER_BIN}" ]]; then
  echo "[bulkde-e2e] ERROR: rs-reborn runner not found on PATH (needs: rvx run)." >&2
  echo "[bulkde-e2e] Hint: go install github.com/rainoffallingstar/rs-reborn/cmd/rvx@latest" >&2
  exit 2
fi

runner_help="$("${RUNNER_BIN}" --help 2>&1 || true)"
if ! echo "${runner_help}" | grep -Eqi "(rs|rvx) run \\[flags\\].*script\\.r"; then
  runner_help="$("${RUNNER_BIN}" run --help 2>&1 || true)"
fi
if ! echo "${runner_help}" | grep -Eqi "(rs|rvx) run \\[flags\\].*script\\.r"; then
  echo "[bulkde-e2e] ERROR: runner found but does not look like rs-reborn (missing: rvx run [flags] path/to/script.R ...)." >&2
  echo "[bulkde-e2e] Hint: go install github.com/rainoffallingstar/rs-reborn/cmd/rvx@latest" >&2
  exit 2
fi

RAW_GZ="${WORK_DIR}/${GSE_ID}.counts.txt.gz"
COUNTS_TSV="${WORK_DIR}/${GSE_ID}.counts.subset.tsv"
META_TSV="${WORK_DIR}/${GSE_ID}.meta.tsv"

if [[ ! -s "${RAW_GZ}" ]]; then
  echo "[bulkde-e2e] downloading ${GSE_ID} supplementary counts..." >&2
  curl -sL --fail "${URL}" -o "${RAW_GZ}"
fi

if [[ ! -s "${COUNTS_TSV}" ]]; then
  echo "[bulkde-e2e] preparing subset counts (first ${N_GENES} genes)..." >&2
  python3 - "${RAW_GZ}" "${COUNTS_TSV}" "${N_GENES}" <<'PY'
import gzip, sys

gz_path, out_path, n_genes = sys.argv[1], sys.argv[2], int(sys.argv[3])

with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f, open(out_path, "w", encoding="utf-8") as w:
    header = None
    # Skip leading comment lines, keep the first real header line.
    for line in f:
        if line.startswith("#"):
            continue
        header = line
        break
    if header is None:
        raise SystemExit("no header line found")
    w.write(header)
    for i, line in enumerate(f):
        if i >= n_genes:
            break
        if line.startswith("#"):
            continue
        w.write(line)
PY
fi

if [[ ! -s "${META_TSV}" ]]; then
  echo "[bulkde-e2e] generating meta (random split A/B + batch)..." >&2
  python3 - "${COUNTS_TSV}" "${META_TSV}" <<'PY'
import sys, csv

counts_path, meta_path = sys.argv[1], sys.argv[2]

with open(counts_path, newline="") as f:
    header = f.readline().rstrip("\n").split("\t")
if len(header) < 10:
    raise SystemExit("counts header too short")

gene_id = header[0]
annot = {"Chr","Start","End","Strand","Length"}
sample_cols = [c for c in header[1:] if c not in annot]
if len(sample_cols) < 4:
    raise SystemExit("not enough sample columns detected")

with open(meta_path, "w", newline="") as g:
    w = csv.writer(g, delimiter="\t")
    w.writerow(["sample","group","batch"])
    half = len(sample_cols)//2
    for i, s in enumerate(sample_cols):
        grp = "A" if i < half else "B"
        batch = "b1" if i % 2 == 0 else "b2"
        w.writerow([s, grp, batch])
PY
fi

echo "[bulkde-e2e] building bulkde..." >&2
cd "${ROOT_DIR}"
GOCACHE=/tmp/go-build go build -o bulkde

echo "[bulkde-e2e] running bulkde..." >&2
args=(
  run
  --counts "${COUNTS_TSV}"
  --out "${OUT_DIR}"
  --group-from meta
  --meta "${META_TSV}"
  --case A
  --control B
  --covariates batch
  --methods "${METHODS}"
  --gene-id-col Geneid
  --annot-cols Chr,Start,End,Strand,Length
  --cache-dir "${CACHE_DIR}"
)
if [[ "${NO_INSTALL}" == "1" ]]; then
  args+=( --no-install )
fi

./bulkde "${args[@]}"

echo "[bulkde-e2e] verifying outputs..." >&2
req=(
  sample_meta.tsv
  run_config.tsv
  design_matrix.tsv
  contrast.tsv
  gene_count.filtered.tsv
  filter_summary.tsv
  de_summary.tsv
  session_info.txt
)
for f in "${req[@]}"; do
  test -s "${OUT_DIR}/${f}"
done

echo "[bulkde-e2e] OK: ${GSE_ID} end-to-end passed. out=${OUT_DIR}" >&2
