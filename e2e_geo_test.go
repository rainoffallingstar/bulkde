package main

import (
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

func TestE2EGEO_GSE116899(t *testing.T) {
	if os.Getenv("BULKDE_E2E") != "1" {
		t.Skip("set BULKDE_E2E=1 to run GEO end-to-end test (downloads data)")
	}

	wd, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	script := filepath.Join(wd, "scripts", "e2e_geo_gse116899.sh")
	cmd := exec.Command("bash", script)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	cmd.Env = os.Environ()
	if err := cmd.Run(); err != nil {
		t.Fatalf("e2e failed: %v", err)
	}
}

