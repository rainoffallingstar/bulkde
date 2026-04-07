package main

import (
	"crypto/sha256"
	_ "embed"
	"encoding/hex"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/rainoffallingstar/rs-reborn/pkg/runner"
)

//go:embed scripts/bulkde.R
var rScript string

type config struct {
	countsFile      string
	outDir          string
	geneIDCol       string
	annotPrefix     string
	annotCols       string
	geneNameCol     string

	groupFrom      string
	metaFile       string
	metaSampleCol  string
	groupCol       string
	caseLabel      string
	controlLabel   string
	covariates     string

	markerGene      string
	markerField     string
	markerThreshold string
	markerOp        string

	cacheDir        string
	compatRLib      string
	noInstall       bool

	minCount        string
	minSamples      string
	varDropQuantile string

	methods      string
	fdrThreshold string
	lfcThreshold string

	// Compat flags
	compatCount string
}

func (c *config) validate() error {
	if c.countsFile == "" {
		c.countsFile = c.compatCount
	}
	if c.countsFile == "" {
		return fmt.Errorf("--counts/--count is required")
	}
	if c.outDir == "" {
		return fmt.Errorf("--out is required")
	}
	if c.markerGene == "" {
		c.markerGene = "FOLH1"
	}
	if c.groupFrom == "" {
		c.groupFrom = "auto"
	}
	if c.metaSampleCol == "" {
		c.metaSampleCol = "sample"
	}
	if c.groupCol == "" {
		c.groupCol = "group"
	}
	if c.geneIDCol == "" {
		c.geneIDCol = "ID"
	}
	if c.annotPrefix == "" {
		c.annotPrefix = "gene_"
	}
	if c.geneNameCol == "" {
		c.geneNameCol = "gene_name"
	}
	if c.markerField == "" {
		c.markerField = "gene_name"
	}
	if c.markerThreshold == "" {
		c.markerThreshold = "0"
	}
	if c.markerOp == "" {
		c.markerOp = "gt"
	}
	if c.minCount == "" {
		c.minCount = "10"
	}
	if c.minSamples == "" {
		c.minSamples = "3"
	}
	if c.varDropQuantile == "" {
		c.varDropQuantile = "0.25"
	}
	if c.methods == "" {
		c.methods = "limma,deseq2,edger"
	}
	if c.fdrThreshold == "" {
		c.fdrThreshold = "0.05"
	}
	if c.lfcThreshold == "" {
		c.lfcThreshold = "1"
	}
	return nil
}

func main() {
	// Support optional subcommand: `bulkde run ...`
	args := os.Args[1:]
	if len(args) > 0 && args[0] == "run" {
		args = args[1:]
	}

	var cfg config
	fs := flag.NewFlagSet("bulkde", flag.ExitOnError)

	// Inputs
	fs.StringVar(&cfg.countsFile, "counts", "", "Counts matrix TSV (gene_id + samples; optional annotation columns)")
	fs.StringVar(&cfg.compatCount, "count", "", "Alias of --counts (compat)")
	fs.StringVar(&cfg.outDir, "out", "", "Output directory")
	fs.StringVar(&cfg.geneIDCol, "gene-id-col", "ID", "Gene ID column name (fallback to first column if missing)")
	fs.StringVar(&cfg.annotPrefix, "annot-prefix", "gene_", "Annotation column prefix (excluded from sample columns)")
	fs.StringVar(&cfg.annotCols, "annot-cols", "", "Comma-separated annotation column names to exclude from sample columns (case-insensitive; e.g. Chr,Start,End,Strand,Length)")
	fs.StringVar(&cfg.geneNameCol, "gene-name-col", "gene_name", "Gene name column (optional; used for marker lookup)")

	// Grouping
	fs.StringVar(&cfg.groupFrom, "group-from", "auto", "Grouping source: auto|meta|marker (auto: meta if --meta provided, else marker)")
	fs.StringVar(&cfg.metaFile, "meta", "", "Sample metadata file (tsv/csv) with at least sample+group")
	fs.StringVar(&cfg.metaSampleCol, "meta-sample-col", "sample", "Sample column name in meta")
	fs.StringVar(&cfg.groupCol, "group-col", "group", "Group column name in meta")
	fs.StringVar(&cfg.caseLabel, "case", "", "Case label in meta group column (case - control)")
	fs.StringVar(&cfg.controlLabel, "control", "", "Control label in meta group column (case - control)")
	fs.StringVar(&cfg.covariates, "covariates", "", "Comma-separated covariate column names from meta (e.g. batch,sex)")

	// Marker grouping options
	fs.StringVar(&cfg.markerGene, "marker", "FOLH1", "Marker gene identifier for marker grouping")
	fs.StringVar(&cfg.markerField, "marker-field", "gene_name", "Marker lookup field: gene_name|gene_id")
	fs.StringVar(&cfg.markerThreshold, "marker-threshold", "0", "Marker threshold for grouping")
	fs.StringVar(&cfg.markerOp, "marker-op", "gt", "Marker op: gt|ge|eq|ne|lt|le (gt means counts > threshold => POS)")

	// Filters
	fs.StringVar(&cfg.minCount, "min-count", "10", "Expression filter: counts >= min-count in >= min-samples samples")
	fs.StringVar(&cfg.minSamples, "min-samples", "3", "Expression filter: counts >= min-count in >= min-samples samples")
	fs.StringVar(&cfg.varDropQuantile, "var-drop-quantile", "0.25", "Variance filter drop quantile (0 disables)")

	// Methods and thresholds
	fs.StringVar(&cfg.methods, "methods", "limma,deseq2,edger", "Methods: limma,deseq2,edger (comma-separated subset allowed)")
	fs.StringVar(&cfg.fdrThreshold, "fdr", "0.05", "Significance threshold for *_sig.tsv (FDR < fdr)")
	fs.StringVar(&cfg.lfcThreshold, "lfc", "1", "Significance threshold for *_sig.tsv (|log2FC| >= lfc)")

	// R execution (rs-reborn)
	fs.StringVar(&cfg.cacheDir, "cache-dir", "", "rs-reborn cache dir (default: <counts_dir>/r_libs)")
	fs.StringVar(&cfg.compatRLib, "r-lib", "", "Alias of --cache-dir (compat)")
	fs.BoolVar(&cfg.noInstall, "no-install", false, "Do not attempt to install R packages (fail if missing)")

	_ = fs.Parse(args)
	if err := cfg.validate(); err != nil {
		fmt.Fprintln(os.Stderr, "ERROR:", err)
		fs.Usage()
		os.Exit(2)
	}

	// Resolve paths early for nicer logging/errors.
	countAbs, err := filepath.Abs(cfg.countsFile)
	if err != nil {
		fmt.Fprintln(os.Stderr, "ERROR: resolve --counts:", err)
		os.Exit(2)
	}
	outAbs, err := filepath.Abs(cfg.outDir)
	if err != nil {
		fmt.Fprintln(os.Stderr, "ERROR: resolve --out:", err)
		os.Exit(2)
	}
	if err := os.MkdirAll(outAbs, 0o755); err != nil {
		fmt.Fprintln(os.Stderr, "ERROR: create output dir:", err)
		os.Exit(2)
	}

	cacheDir := cfg.cacheDir
	if cacheDir == "" {
		cacheDir = cfg.compatRLib
	}
	if cacheDir == "" {
		cacheDir = filepath.Join(filepath.Dir(countAbs), "r_libs")
	}
	cacheDirAbs, err := filepath.Abs(cacheDir)
	if err != nil {
		fmt.Fprintln(os.Stderr, "ERROR: resolve --cache-dir:", err)
		os.Exit(2)
	}
	if err := os.MkdirAll(cacheDirAbs, 0o755); err != nil {
		fmt.Fprintln(os.Stderr, "ERROR: create cache dir:", err)
		os.Exit(2)
	}

	scriptPath, err := materializeEmbeddedScript("bulkde.R", rScript)
	if err != nil {
		fmt.Fprintln(os.Stderr, "ERROR: materialize embedded R script:", err)
		os.Exit(2)
	}

	fmt.Fprintf(os.Stderr, "bulkde: counts=%s out=%s group-from=%s methods=%s\n", countAbs, outAbs, cfg.groupFrom, cfg.methods)

	// Pass config via env (keeps R argv stable).
	env := map[string]string{
		"BULKDE_GENE_ID_COL":        cfg.geneIDCol,
		"BULKDE_ANNOT_PREFIX":       cfg.annotPrefix,
		"BULKDE_ANNOT_COLS":         cfg.annotCols,
		"BULKDE_GENE_NAME_COL":      cfg.geneNameCol,
		"BULKDE_GROUP_FROM":         cfg.groupFrom,
		"BULKDE_META":               cfg.metaFile,
		"BULKDE_META_SAMPLE_COL":    cfg.metaSampleCol,
		"BULKDE_GROUP_COL":          cfg.groupCol,
		"BULKDE_CASE":               cfg.caseLabel,
		"BULKDE_CONTROL":            cfg.controlLabel,
		"BULKDE_COVARIATES":         cfg.covariates,
		"BULKDE_MARKER":             cfg.markerGene,
		"BULKDE_MARKER_FIELD":       cfg.markerField,
		"BULKDE_MARKER_THRESHOLD":   cfg.markerThreshold,
		"BULKDE_MARKER_OP":          cfg.markerOp,
		"BULKDE_MIN_COUNT":          cfg.minCount,
		"BULKDE_MIN_SAMPLES":        cfg.minSamples,
		"BULKDE_VAR_DROP_QUANTILE":  cfg.varDropQuantile,
		"BULKDE_METHODS":            cfg.methods,
		"BULKDE_FDR":                cfg.fdrThreshold,
		"BULKDE_LFC":                cfg.lfcThreshold,
	}

	exclude := excludeDepsForMethods(cfg.methods)

	if err := withEnv(env, func() error {
		return runner.Run(runner.RunOptions{
			ScriptPath:    scriptPath,
			ScriptArgs:    []string{countAbs, outAbs},
			CacheDir:      cacheDirAbs,
			SkipInstall:   cfg.noInstall,
			ExcludeDeps:   exclude,
			Stdout:        os.Stdout,
			Stderr:        os.Stderr,
			AutoInstallR:  false,
			Verbose:       false,
			Locked:        false,
			Frozen:        false,
			BootstrapToolchain: false,
		})
	}); err != nil {
		fmt.Fprintln(os.Stderr, "ERROR: differential analysis failed:", err)
		os.Exit(1)
	}
}

func withEnv(kv map[string]string, fn func() error) error {
	// Save old values and restore after run (runner inherits from os.Environ()).
	old := make(map[string]*string, len(kv))
	for k, v := range kv {
		if v == "" {
			continue
		}
		if ov, ok := os.LookupEnv(k); ok {
			tmp := ov
			old[k] = &tmp
		} else {
			old[k] = nil
		}
		_ = os.Setenv(k, v)
	}
	defer func() {
		for k, ov := range old {
			if ov == nil {
				_ = os.Unsetenv(k)
			} else {
				_ = os.Setenv(k, *ov)
			}
		}
	}()
	return fn()
}

func materializeEmbeddedScript(filename string, content string) (string, error) {
	sum := sha256.Sum256([]byte(content))
	hash := hex.EncodeToString(sum[:])

	base := os.TempDir()
	if d, err := os.UserCacheDir(); err == nil && d != "" {
		base = d
	}

	dir := filepath.Join(base, "bulkde", "scripts", hash)
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return "", err
	}
	path := filepath.Join(dir, filename)

	// If exists and non-empty, reuse.
	if st, err := os.Stat(path); err == nil && st.Size() > 0 {
		return path, nil
	}

	// Atomic-ish write.
	tmp := filepath.Join(dir, filename+".tmp")
	if err := os.WriteFile(tmp, []byte(content), 0o755); err != nil {
		return "", err
	}
	if err := os.Rename(tmp, path); err != nil {
		// If rename fails (e.g. cross-device), fallback to direct write.
		_ = os.Remove(tmp)
		if err2 := os.WriteFile(path, []byte(content), 0o755); err2 != nil {
			return "", err2
		}
	}
	return path, nil
}

func excludeDepsForMethods(methodsRaw string) []string {
	m := parseMethods(methodsRaw)
	keepLimma := m["limma"]
	keepDESeq2 := m["deseq2"]
	// edgeR is always needed for CPM filtering.
	var exclude []string
	if !keepLimma {
		exclude = append(exclude, "limma")
	}
	if !keepDESeq2 {
		exclude = append(exclude, "DESeq2")
	}
	return exclude
}

func parseMethods(methodsRaw string) map[string]bool {
	out := map[string]bool{}
	for _, p := range strings.Split(methodsRaw, ",") {
		p = strings.TrimSpace(strings.ToLower(p))
		if p == "" {
			continue
		}
		out[p] = true
	}
	return out
}
