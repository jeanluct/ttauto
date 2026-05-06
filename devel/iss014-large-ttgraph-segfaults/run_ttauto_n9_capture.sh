#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

TTAUTO_BIN="${TTAUTO_BIN:-${REPO_DIR}/examples/ttauto}"
INPUT_FILE="${INPUT_FILE:-${SCRIPT_DIR}/ttauto_n9.in}"
CAPTURE_MODE="${CAPTURE_MODE:-fast}"   # fast | gdb | rr
RUN_ROOT="${RUN_ROOT:-${SCRIPT_DIR}/crash_runs}"
MAX_SECONDS="${MAX_SECONDS:-0}"         # 0 means no timeout
ADD_VERTEX_PROGRESS_EVERY="${ADD_VERTEX_PROGRESS_EVERY:-0}"
ADD_VERTEX_PROGRESS_CALLS_EVERY="${ADD_VERTEX_PROGRESS_CALLS_EVERY:-0}"

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
Usage: devel/iss014-large-ttgraph-segfaults/run_ttauto_n9_capture.sh

Environment knobs:
  TTAUTO_BIN   Path to ttauto executable (default: examples/ttauto)
  INPUT_FILE   Path to replay input file (default: devel/iss014-large-ttgraph-segfaults/ttauto_n9.in)
  CAPTURE_MODE fast | gdb | rr (default: fast)
  RUN_ROOT     Output directory root (default: devel/iss014-large-ttgraph-segfaults/crash_runs)
  MAX_SECONDS  Stop run after N seconds (0: no timeout)
  ADD_VERTEX_PROGRESS_EVERY       Emit progress every N vertices (0: off)
  ADD_VERTEX_PROGRESS_CALLS_EVERY Emit progress every N add_vertex calls (0: off)

Examples:
  CAPTURE_MODE=fast devel/iss014-large-ttgraph-segfaults/run_ttauto_n9_capture.sh
  CAPTURE_MODE=gdb  devel/iss014-large-ttgraph-segfaults/run_ttauto_n9_capture.sh
EOF
  exit 0
fi

if ! [[ "${MAX_SECONDS}" =~ ^[0-9]+$ ]]; then
  printf "Error: MAX_SECONDS must be a nonnegative integer (got %s).\n" "${MAX_SECONDS}" >&2
  exit 1
fi

if ! [[ "${ADD_VERTEX_PROGRESS_EVERY}" =~ ^[0-9]+$ ]]; then
  printf "Error: ADD_VERTEX_PROGRESS_EVERY must be a nonnegative integer (got %s).\n" \
    "${ADD_VERTEX_PROGRESS_EVERY}" >&2
  exit 1
fi

if ! [[ "${ADD_VERTEX_PROGRESS_CALLS_EVERY}" =~ ^[0-9]+$ ]]; then
  printf "Error: ADD_VERTEX_PROGRESS_CALLS_EVERY must be a nonnegative integer (got %s).\n" \
    "${ADD_VERTEX_PROGRESS_CALLS_EVERY}" >&2
  exit 1
fi

if [[ ! -x "${TTAUTO_BIN}" ]]; then
  printf "Error: executable not found: %s\n" "${TTAUTO_BIN}" >&2
  printf "Build first, e.g. cmake -S . -B build && cmake --build build -j\n" >&2
  exit 1
fi

if [[ ! -f "${INPUT_FILE}" ]]; then
  printf "Error: input file not found: %s\n" "${INPUT_FILE}" >&2
  exit 1
fi

mkdir -p "${RUN_ROOT}"
stamp="$(date +%Y%m%d_%H%M%S)"
run_dir="${RUN_ROOT}/${stamp}_${CAPTURE_MODE}"
mkdir -p "${run_dir}"

run_log="${run_dir}/run.log"
meta_log="${run_dir}/meta.txt"
gdb_log="${run_dir}/gdb.log"
core_log="${run_dir}/core_files.txt"

copy_cores_from_tmp() {
  local copied=0
  shopt -s nullglob
  for f in /tmp/core /tmp/core.*; do
    [[ -f "${f}" ]] || continue
    cp -n "${f}" "${run_dir}/" || true
    copied=1
  done
  shopt -u nullglob
  return "${copied}"
}

ulimit -c unlimited || true

{
  printf "timestamp=%s\n" "$(date -Is)"
  printf "mode=%s\n" "${CAPTURE_MODE}"
  printf "repo=%s\n" "${REPO_DIR}"
  printf "bin=%s\n" "${TTAUTO_BIN}"
  printf "input=%s\n" "${INPUT_FILE}"
  printf "cwd=%s\n" "${run_dir}"
  git -C "${REPO_DIR}" rev-parse HEAD 2>/dev/null | sed 's/^/git_head=/'
} > "${meta_log}"

printf "Capture directory: %s\n" "${run_dir}"

status=0
case "${CAPTURE_MODE}" in
  fast)
    # Minimal-overhead capture: deterministic replay + full stdout/stderr + core dump.
    (
      cd "${run_dir}"
      if (( MAX_SECONDS > 0 )); then
        timeout --signal=INT --kill-after=10s "${MAX_SECONDS}" env \
          TTAUTO_ADD_VERTEX_PROGRESS_EVERY="${ADD_VERTEX_PROGRESS_EVERY}" \
          TTAUTO_ADD_VERTEX_PROGRESS_CALLS_EVERY="${ADD_VERTEX_PROGRESS_CALLS_EVERY}" \
          "${TTAUTO_BIN}" < "${INPUT_FILE}" > "${run_log}" 2>&1
      else
        env \
          TTAUTO_ADD_VERTEX_PROGRESS_EVERY="${ADD_VERTEX_PROGRESS_EVERY}" \
          TTAUTO_ADD_VERTEX_PROGRESS_CALLS_EVERY="${ADD_VERTEX_PROGRESS_CALLS_EVERY}" \
          "${TTAUTO_BIN}" < "${INPUT_FILE}" > "${run_log}" 2>&1
      fi
    ) || status=$?
    ;;

  gdb)
    # Heavier capture: gdb batch run with automatic full backtrace.
    (
      cd "${run_dir}"
      gdb_cmd=(
        gdb -batch
        -ex "set pagination off"
      )
      if (( MAX_SECONDS > 0 )); then
        gdb_cmd+=(
          -ex "python import threading,gdb; threading.Timer(${MAX_SECONDS}, lambda: gdb.execute('interrupt')).start()"
        )
      fi
      gdb_cmd+=(
        -ex run
        -ex "thread apply all bt full"
        -ex "info threads"
        --args env
        TTAUTO_ADD_VERTEX_PROGRESS_EVERY="${ADD_VERTEX_PROGRESS_EVERY}"
        TTAUTO_ADD_VERTEX_PROGRESS_CALLS_EVERY="${ADD_VERTEX_PROGRESS_CALLS_EVERY}"
        "${TTAUTO_BIN}"
      )
      "${gdb_cmd[@]}" < "${INPUT_FILE}" > "${gdb_log}" 2>&1
    ) || status=$?

    # gdb often exits 0 even on inferior crash/interrupt. Derive status from log.
    if grep -q "Program received signal SIGSEGV" "${gdb_log}"; then
      status=139
    elif (( MAX_SECONDS > 0 )) && grep -q "Program received signal SIGINT" "${gdb_log}"; then
      status=124
    fi
    ;;

  rr)
    if ! command -v rr >/dev/null 2>&1; then
      printf "Error: CAPTURE_MODE=rr but rr is not installed.\n" >&2
      exit 1
    fi
    (
      cd "${run_dir}"
      if (( MAX_SECONDS > 0 )); then
        timeout --signal=INT --kill-after=10s "${MAX_SECONDS}" \
          rr record env \
          TTAUTO_ADD_VERTEX_PROGRESS_EVERY="${ADD_VERTEX_PROGRESS_EVERY}" \
          TTAUTO_ADD_VERTEX_PROGRESS_CALLS_EVERY="${ADD_VERTEX_PROGRESS_CALLS_EVERY}" \
          "${TTAUTO_BIN}" < "${INPUT_FILE}" > "${run_log}" 2>&1
      else
        rr record env \
          TTAUTO_ADD_VERTEX_PROGRESS_EVERY="${ADD_VERTEX_PROGRESS_EVERY}" \
          TTAUTO_ADD_VERTEX_PROGRESS_CALLS_EVERY="${ADD_VERTEX_PROGRESS_CALLS_EVERY}" \
          "${TTAUTO_BIN}" < "${INPUT_FILE}" > "${run_log}" 2>&1
      fi
    ) || status=$?
    ;;

  *)
    printf "Error: unsupported CAPTURE_MODE=%s (use fast|gdb|rr).\n" "${CAPTURE_MODE}" >&2
    exit 1
    ;;
esac

{
  printf "exit_status=%s\n" "${status}"
  printf "core_ulimit=%s\n" "$(ulimit -c 2>/dev/null || printf unknown)"
  printf "max_seconds=%s\n" "${MAX_SECONDS}"
  printf "progress_vertices_every=%s\n" "${ADD_VERTEX_PROGRESS_EVERY}"
  printf "progress_calls_every=%s\n" "${ADD_VERTEX_PROGRESS_CALLS_EVERY}"
} >> "${meta_log}"

(
  cd "${run_dir}"
  found=0
  for pat in core core.*; do
    for f in $pat; do
      [[ -e "${f}" ]] || continue
      printf "%s\n" "${f}"
      found=1
    done
  done
  if [[ "${found}" -eq 0 ]]; then
    printf "(no core files found)\n"
  fi
) > "${core_log}"

if [[ "${status}" -ne 0 ]]; then
  if [[ "${status}" -eq 124 ]]; then
    printf "Run reached timeout after %ss.\n" "${MAX_SECONDS}"
  else
    printf "Run exited with status %s (expected for crash reproduction).\n" "${status}"
  fi
  copy_cores_from_tmp || true
else
  printf "Run completed without nonzero exit.\n"
fi

printf "Artifacts:\n"
printf "  %s\n" "${meta_log}"
if [[ "${CAPTURE_MODE}" == "gdb" ]]; then
  printf "  %s\n" "${gdb_log}"
else
  printf "  %s\n" "${run_log}"
fi
printf "  %s\n" "${core_log}"

exit "${status}"
