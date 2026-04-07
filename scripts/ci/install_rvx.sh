#!/usr/bin/env bash
set -euo pipefail

# Install the latest available rs-reborn runner binary (rvx) for the current platform.
# Preference order:
#  1) Download from GitHub Releases "latest" assets (fast, stable)
#  2) Fallback: go install from source (@latest)
#
# Output:
#  - Ensures `rvx` is on PATH (caller should export PATH if needed)

REPO="${RVX_REPO:-rainoffallingstar/rs-reborn}"
OUT_DIR="${RVX_INSTALL_DIR:-${PWD}/.ci-bin}"
FALLBACK_GO_INSTALL="${RVX_FALLBACK_GO_INSTALL:-1}"

mkdir -p "${OUT_DIR}"

if command -v rvx >/dev/null 2>&1; then
  echo "[install_rvx] rvx already present: $(command -v rvx)" >&2
  exit 0
fi

os="$(uname -s | tr '[:upper:]' '[:lower:]')"
arch="$(uname -m)"
case "${arch}" in
  x86_64|amd64) arch="amd64" ;;
  aarch64|arm64) arch="arm64" ;;
esac

api="https://api.github.com/repos/${REPO}/releases/latest"
tmp_json="$(mktemp)"
tmp_dir="$(mktemp -d)"
cleanup() { rm -f "${tmp_json}"; rm -rf "${tmp_dir}"; }
trap cleanup EXIT

curl_args=(-fsSL -H "Accept: application/vnd.github+json")
if [[ -n "${GITHUB_TOKEN:-}" ]]; then
  curl_args+=(-H "Authorization: Bearer ${GITHUB_TOKEN}")
fi

echo "[install_rvx] probing latest release: ${api}" >&2
if curl "${curl_args[@]}" "${api}" -o "${tmp_json}"; then
  if command -v jq >/dev/null 2>&1; then
    # Try common asset naming patterns.
    # We accept either an archive that contains `rvx` or a raw binary asset.
    name_regex="rvx.*${os}.*${arch}"
    asset_url="$(jq -r --arg re "${name_regex}" '
      (.assets // [])
      | map(select(.name|test($re; "i")))
      | map(select(.browser_download_url != null))
      | (map(select(.name|test("\\.(tar\\.gz|tgz|zip)$"; "i"))) + map(select(.name|test("\\.(tar\\.gz|tgz|zip)$"; "i")|not)))
      | .[0].browser_download_url // empty
    ' "${tmp_json}")"

    if [[ -n "${asset_url}" ]]; then
      echo "[install_rvx] downloading: ${asset_url}" >&2
      asset_path="${tmp_dir}/asset"
      curl -fsSL "${asset_url}" -o "${asset_path}"

      # Determine type by filename if possible; otherwise by magic.
      asset_name="$(basename "${asset_url}")"
      if [[ "${asset_name}" =~ \.zip$ ]]; then
        unzip -q "${asset_path}" -d "${tmp_dir}/unzip"
        found="$(find "${tmp_dir}/unzip" -type f -name rvx -o -name rvx.exe | head -n 1 || true)"
        if [[ -z "${found}" ]]; then
          echo "[install_rvx] zip downloaded but rvx not found inside" >&2
        else
          cp "${found}" "${OUT_DIR}/rvx"
          chmod +x "${OUT_DIR}/rvx"
          echo "[install_rvx] installed to ${OUT_DIR}/rvx" >&2
          exit 0
        fi
      elif [[ "${asset_name}" =~ \.(tar\.gz|tgz)$ ]]; then
        tar -xzf "${asset_path}" -C "${tmp_dir}/untar"
        found="$(find "${tmp_dir}/untar" -type f -name rvx -o -name rvx.exe | head -n 1 || true)"
        if [[ -z "${found}" ]]; then
          echo "[install_rvx] tarball downloaded but rvx not found inside" >&2
        else
          cp "${found}" "${OUT_DIR}/rvx"
          chmod +x "${OUT_DIR}/rvx"
          echo "[install_rvx] installed to ${OUT_DIR}/rvx" >&2
          exit 0
        fi
      else
        # Assume it's a raw binary.
        cp "${asset_path}" "${OUT_DIR}/rvx"
        chmod +x "${OUT_DIR}/rvx"
        echo "[install_rvx] installed raw asset to ${OUT_DIR}/rvx" >&2
        exit 0
      fi
    else
      echo "[install_rvx] no matching release asset for ${os}/${arch}; falling back" >&2
    fi
  else
    echo "[install_rvx] jq not available; skipping release-asset install" >&2
  fi
else
  echo "[install_rvx] GitHub API request failed; falling back" >&2
fi

if [[ "${FALLBACK_GO_INSTALL}" == "1" ]]; then
  echo "[install_rvx] go install fallback: github.com/${REPO}/cmd/rvx@latest" >&2
  GOBIN="${OUT_DIR}" go install "github.com/${REPO}/cmd/rvx@latest"
  test -x "${OUT_DIR}/rvx"
  echo "[install_rvx] installed via go to ${OUT_DIR}/rvx" >&2
  exit 0
fi

echo "[install_rvx] ERROR: could not install rvx (no matching release asset; go fallback disabled)" >&2
exit 2

