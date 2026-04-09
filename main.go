package main

import (
	"crypto/sha256"
	_ "embed"
	"encoding/hex"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
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
	rvxPath         string
	compatRsPath    string

	minCount        string
	minSamples      string
	varDropQuantile string

	methods      string
	fdrThreshold string
	lfcThreshold string
	chipMode     bool

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
	if c.chipMode {
		c.methods = "limma"
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
	fs.BoolVar(&cfg.chipMode, "chip-mode", false, "Treat input as preprocessed microarray/chip expression matrix; forces limma-only and skips RNA-seq-specific filtering")

	// R execution (rs-reborn)
	fs.StringVar(&cfg.cacheDir, "cache-dir", "", "rs-reborn cache dir (default: <counts_dir>/r_libs)")
	fs.StringVar(&cfg.compatRLib, "r-lib", "", "Alias of --cache-dir (compat)")
	fs.BoolVar(&cfg.noInstall, "no-install", false, "Do not attempt to install R packages (fail if missing)")
	fs.StringVar(&cfg.rvxPath, "rvx-path", "", "Path to rs-reborn CLI binary (rvx; default: find `rvx` from PATH)")
	fs.StringVar(&cfg.compatRsPath, "rs-path", "", "Alias of --rvx-path (compat; old binary name)")

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
		"BULKDE_CHIP_MODE":          bool01(cfg.chipMode),
		"BULKDE_NO_INSTALL":         bool01(cfg.noInstall),
	}

	exclude := excludeDepsForMethods(cfg.methods, cfg.chipMode)

	runnerPath := cfg.rvxPath
	if runnerPath == "" {
		runnerPath = cfg.compatRsPath
	}
	if err := runWithRebornCLI(runWithRebornCLIOptions{
		RunnerPath: runnerPath,
		CacheDir:   cacheDirAbs,
		SkipInstall: cfg.noInstall,
		Exclude:    exclude,
		ScriptPath: scriptPath,
		ScriptArgs: []string{countAbs, outAbs},
		Env:        env,
	}); err != nil {
		fmt.Fprintln(os.Stderr, "ERROR: differential analysis failed:", err)
		os.Exit(1)
	}
}

type runWithRebornCLIOptions struct {
	RunnerPath  string
	CacheDir    string
	SkipInstall bool
	Exclude     []string
	ScriptPath  string
	ScriptArgs  []string
	Env         map[string]string
}

func runWithRebornCLI(opt runWithRebornCLIOptions) error {
	runnerPath := strings.TrimSpace(opt.RunnerPath)
	runnerName := ""
	if runnerPath == "" {
		if p, err := exec.LookPath("rvx"); err == nil {
			runnerPath = p
			runnerName = "rvx"
		} else if p, err := exec.LookPath("rs"); err == nil {
			runnerPath = p
			runnerName = "rs"
		} else {
			return fmt.Errorf("rvx CLI not found in PATH; please install rs-reborn runner (e.g. `go install github.com/rainoffallingstar/rs-reborn/cmd/rvx@latest`) or pass --rvx-path")
		}
	}
	if runnerName == "" {
		runnerName = filepath.Base(runnerPath)
	}
	if ok, why := looksLikeRebornRunner(runnerPath); !ok {
		return fmt.Errorf("found runner at %s but it does not look like rs-reborn (missing `%s run` usage): %s. Please install rs-reborn runner and ensure it is first on PATH, or pass --rvx-path to the correct binary", runnerPath, runnerName, why)
	}

	args := []string{"run", "--cache-dir", opt.CacheDir}
	if opt.SkipInstall {
		args = append(args, "--no-install")
	}
	for _, ex := range opt.Exclude {
		if strings.TrimSpace(ex) == "" {
			continue
		}
		args = append(args, "--exclude", ex)
	}
	args = append(args, opt.ScriptPath)
	args = append(args, opt.ScriptArgs...)

	fmt.Fprintf(os.Stderr, "bulkde: exec: %s %s\n", runnerPath, strings.Join(args, " "))

	cmd := exec.Command(runnerPath, args...)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	// Inherit environment, then override with BULKDE_* variables.
	env := os.Environ()
	for k, v := range opt.Env {
		if strings.TrimSpace(k) == "" || strings.TrimSpace(v) == "" {
			continue
		}
		env = append(env, k+"="+v)
	}
	cmd.Env = env

	if err := cmd.Run(); err != nil {
		return err
	}
	return nil
}

func looksLikeRebornRunner(runnerPath string) (bool, string) {
	// Note: many systems already ship a different `/usr/bin/rs` (column formatting tool),
	// so we validate presence of the "run" subcommand and the script-based usage shape.
	//
	// rs-reborn typically contains a line like:
	//   rvx run [flags] path/to/script.R [script args...]
	// or:
	//   rs run [flags] path/to/script.R [script args...]
	// Some versions print this under `rs --help` (top-level), some under `rs run --help`,
	// and some exit non-zero for `--help`. We accept any of these as long as the text matches.
	matches := func(s string) bool {
		h := strings.ToLower(s)
		if strings.Contains(h, "rs manages a lightweight per-script r library") {
			return true
		}
		if strings.Contains(h, " run [flags]") && strings.Contains(h, "script.r") {
			return true
		}
		if strings.Contains(h, "usage: rs run") {
			return true
		}
		if strings.Contains(h, "usage: rvx run") {
			return true
		}
		return false
	}

	try := func(args ...string) (bool, string) {
		cmd := exec.Command(runnerPath, args...)
		out, err := cmd.CombinedOutput()
		txt := string(out)
		if matches(txt) {
			return true, ""
		}
		if err != nil {
			msg := strings.TrimSpace(txt)
			if msg == "" {
				msg = err.Error()
			}
			return false, msg
		}
		return false, strings.TrimSpace(txt)
	}

	if ok, why := try("--help"); ok {
		return true, ""
	} else if why != "" {
		// Continue to next probe; collect last message for debugging.
	}
	if ok, why := try("run", "--help"); ok {
		return true, ""
	} else if why != "" {
		return false, why
	}
	return false, "missing rs-reborn markers in help output (expected `rs run [flags] path/to/script.R`)"
}

func bool01(b bool) string {
	if b {
		return "1"
	}
	return "0"
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

func excludeDepsForMethods(methodsRaw string, chipMode bool) []string {
	m := parseMethods(methodsRaw)
	keepLimma := m["limma"]
	keepDESeq2 := m["deseq2"]
	var exclude []string
	if !keepLimma {
		exclude = append(exclude, "limma")
	}
	if !keepDESeq2 {
		exclude = append(exclude, "DESeq2")
	}
	// edgeR is always needed for the RNA-seq path because filtering and limma-voom
	// both rely on DGEList/cpm utilities. Only chip-mode can safely exclude it.
	if chipMode {
		exclude = append(exclude, "edgeR")
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
