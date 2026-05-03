#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TTAUTO_BIN="${TTAUTO_BIN:-${SCRIPT_DIR}/ttauto}"
MAX_PATH_LEN="${MAX_PATH_LEN:-5}"
OUTPUT_MD="${OUTPUT_MD:-${SCRIPT_DIR}/ttauto_scan_strata.md}"
SHOW_PROGRESS="${SHOW_PROGRESS:-1}"
N_MIN="${N_MIN:-3}"
N_MAX="${N_MAX:-7}"

# Per-(puncture,stratum) max path lengths.
#
# Key format is "n:stratum" (e.g. "6:5").
#
# Defaults are hardwired here for easy manual editing.
# For n=7 all known strata default to 4.
declare -A MAX_PATH_LEN_BY_NS=(
  ["3:1"]=4

  ["4:1"]=4
  ["4:2"]=4

  ["5:1"]=4
  ["5:2"]=4
  ["5:3"]=4
  ["5:4"]=4

  ["6:1"]=4
  ["6:2"]=4
  ["6:3"]=5
  ["6:4"]=4
  ["6:5"]=6
  ["6:6"]=4
  ["6:7"]=8

  ["7:1"]=4
  ["7:2"]=4
  ["7:3"]=4
  ["7:4"]=4
  ["7:5"]=4
  ["7:6"]=5
  ["7:7"]=7
  ["7:8"]=4
  ["7:9"]=4
  ["7:10"]=4
  ["7:11"]=8
  ["7:12"]=8
)

progress() {
  if [[ "${SHOW_PROGRESS}" == "1" ]]; then
    printf "%s\n" "$*" >&2
  fi
}

max_len_for_case() {
  local n="$1"
  local stratum="$2"
  local key="${n}:${stratum}"

  if [[ -v MAX_PATH_LEN_BY_NS["${key}"] ]]; then
    printf "%s\n" "${MAX_PATH_LEN_BY_NS["${key}"]}"
    return
  fi

  printf "%s\n" "${MAX_PATH_LEN}"
}

if [[ ! -x "${TTAUTO_BIN}" ]]; then
  printf "Error: executable not found: %s\n" "${TTAUTO_BIN}" >&2
  printf "Build first, e.g. cmake -S . -B build && cmake --build build -j\n" >&2
  exit 1
fi

if ! [[ "${N_MIN}" =~ ^[0-9]+$ && "${N_MAX}" =~ ^[0-9]+$ ]]; then
  printf "Error: N_MIN and N_MAX must be integers (got N_MIN=%s N_MAX=%s)\n" "${N_MIN}" "${N_MAX}" >&2
  exit 1
fi

if (( N_MIN < 3 || N_MAX < 3 || N_MIN > N_MAX )); then
  printf "Error: require 3 <= N_MIN <= N_MAX (got N_MIN=%s N_MAX=%s)\n" "${N_MIN}" "${N_MAX}" >&2
  exit 1
fi

extract_strata_lines() {
  local n="$1"
  local out

  out="$(printf "%s\n1\nn\n1\n0\n" "${n}" | "${TTAUTO_BIN}" 2>/dev/null || true)"

  printf "%s\n" "${out}" | awk '/^[[:space:]]*[0-9]+\)[[:space:]]/{print}'
}

parse_run_metrics() {
  local out="$1"

  local main_size total_size other_size
  local min_dil shortest_minimiser shortest_overall

  main_size="$(printf "%s\n" "${out}" | awk -F'\t' '
    /^[[:space:]]*1[[:space:]]*\t/ {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2);
      print $2;
      exit;
    }
  ')"

  total_size="$(printf "%s\n" "${out}" | awk -F'\t' '
    /^[[:space:]]*total[[:space:]]*\t/ {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2);
      print $2;
      exit;
    }
  ')"

  if [[ -n "${main_size}" && -n "${total_size}" ]]; then
    other_size="$((total_size - main_size))"
  else
    other_size="NA"
  fi

  read -r min_dil shortest_minimiser shortest_overall <<<"$(printf "%s\n" "${out}" | awk '
    function path_length_from_vertices(inner, cleaned, parts, nfields) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", inner);
      if (inner == "") return 0;

      # Vertices are printed either as space-separated integers (when vertex
      # labels can exceed 9) or as contiguous single digits.
      if (inner ~ /[[:space:]]/) {
        nfields = split(inner, parts, /[[:space:]]+/);
        return (nfields > 0 ? nfields - 1 : 0);
      }

      cleaned = inner;
      gsub(/[[:space:]]/, "", cleaned);
      return (length(cleaned) > 0 ? length(cleaned) - 1 : 0);
    }

    /^[[:space:]]*[0-9]+[[:space:]]+[0-9]+\.[0-9]+/ {
      d = $2;
      if (!have || d < min) {
        have = 1;
        min = d;
        best = 1e9;
      }

      line = $0;
      while (match(line, /\[[^]]+\]/)) {
        bracket = substr(line, RSTART + 1, RLENGTH - 2);
        l = path_length_from_vertices(bracket);

        if (overall == 0 || l < overall) overall = l;
        if (d == min && l < best) best = l;

        line = substr(line, RSTART + RLENGTH);
      }
    }
    END {
      if (!have) {
        print "NA NA NA";
      } else {
        if (best == 1e9) best = 0;
        if (overall == 0) overall = 0;
        printf "%s %d %d\n", min, best, overall;
      }
    }
  ')"

  printf "%s\t%s\t%s\t%s\t%s\n" "${main_size:-NA}" "${other_size}" "${min_dil}" "${shortest_minimiser}" "${shortest_overall}"
}

extract_effective_max_len() {
  local out="$1"
  printf "%s\n" "${out}" | awk '
    /pA candidates found with path length <=/ ||
    /NO pseudo-Anosov candidates found with path length <=/ {
      x = $0;
      gsub(/[^0-9]/, "", x);
      print x;
      exit;
    }
  '
}

run_ttauto_for_stratum() {
  local n="$1"
  local stratum="$2"
  local nstrata="$3"
  local max_len="$4"

  local choose_part=""
  if (( nstrata > 1 )); then
    choose_part="${stratum}\n"
  fi

  # First try assumes "Subgraph to search" is prompted.
  local out
  out="$(printf "%s\n%b" "${n}" "${choose_part}y\nn\n1\n${max_len}\n" | "${TTAUTO_BIN}" 2>/dev/null || true)"

  local used_len
  used_len="$(extract_effective_max_len "${out}")"

  # If the run used a different limit, retry assuming no subgraph prompt
  # (single-subgraph case, where sending "1" would shift inputs).
  if [[ -n "${used_len}" && "${used_len}" != "${max_len}" ]]; then
    out="$(printf "%s\n%b" "${n}" "${choose_part}y\nn\n${max_len}\n" | "${TTAUTO_BIN}" 2>/dev/null || true)"
  fi

  printf "%s\n" "${out}"
}

printf "# ttauto Strata Scan (n=%s..%s)\n\n" "${N_MIN}" "${N_MAX}" > "${OUTPUT_MD}"

for ((n=N_MIN; n<=N_MAX; ++n)); do
  strata_lines="$(extract_strata_lines "${n}")"
  nstrata="$(printf "%s\n" "${strata_lines}" | awk 'NF{c++} END{print c+0}')"

  progress "[n=${n}] found ${nstrata} strata"

  {
    printf "## %s Punctures\n\n" "${n}"
    printf "| Stratum | Singularity Data | Main | non-Main | Min Dil. | Min. Dil Len. | Shortest Len. | Max Len. |\n"
    printf "|---:|---|---:|---:|---:|---:|---:|---:|\n"
  } >> "${OUTPUT_MD}"

  if [[ -z "${strata_lines}" ]]; then
    printf "| NA | NA | NA | NA | NA | NA | NA | NA |\n\n" >> "${OUTPUT_MD}"
    continue
  fi

  while IFS= read -r line; do
    [[ -z "${line}" ]] && continue

    stratum="$(printf "%s\n" "${line}" | sed -E 's/^[[:space:]]*([0-9]+)\).*/\1/')"
    singularity="$(printf "%s\n" "${line}" | sed -E 's/^[[:space:]]*[0-9]+\)[[:space:]]*//')"
    max_len_for_case_val="$(max_len_for_case "${n}" "${stratum}")"

    progress "  -> running stratum ${stratum}: ${singularity} (max path length ${max_len_for_case_val})"

    run_out="$(run_ttauto_for_stratum "${n}" "${stratum}" "${nstrata}" "${max_len_for_case_val}")"

    read -r main_vertices other_vertices_total min_dilatation minimiser_path_len shortest_path_len_overall <<<"$(parse_run_metrics "${run_out}")"

    progress "     main=${main_vertices} other=${other_vertices_total} min_dil=${min_dilatation} min_path=${minimiser_path_len} shortest_any=${shortest_path_len_overall}"

    singularity_md="${singularity//|/\\|}"
    shortest_path_len_display="${shortest_path_len_overall}"
    if [[ "${shortest_path_len_overall}" != "${minimiser_path_len}" ]]; then
      shortest_path_len_display="**${shortest_path_len_overall}**"
    fi

    printf "| %s | %s | %s | %s | %s | %s | %s | %s |\n" \
      "${stratum}" \
      "${singularity_md}" \
      "${main_vertices}" \
      "${other_vertices_total}" \
      "${min_dilatation}" \
      "${minimiser_path_len}" \
      "${shortest_path_len_display}" \
      "${max_len_for_case_val}" >> "${OUTPUT_MD}"
  done <<<"${strata_lines}"

  printf "\n" >> "${OUTPUT_MD}"
done

printf "Wrote markdown report: %s\n" "${OUTPUT_MD}"
