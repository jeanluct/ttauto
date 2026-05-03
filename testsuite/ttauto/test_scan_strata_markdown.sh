#!/usr/bin/env bash

# <LICENSE
#   ttauto: a C++ library for building train track automata
#
#   https://github.com/jeanluct/ttauto
#
#   Copyright (C) 2010-2026  Jean-Luc Thiffeault   <jeanluc@math.wisc.edu>
#                            Erwan Lanneau <erwan.lanneau@ujf-grenoble.fr>
#
#   This file is part of ttauto.
#
#   ttauto is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   ttauto is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ttauto.  If not, see <http://www.gnu.org/licenses/>.
# LICENSE>

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SCAN_SCRIPT="${REPO_DIR}/examples/ttauto_scan_strata.sh"
BASELINE_MD="${REPO_DIR}/examples/ttauto_scan_strata.md"
TTAUTO_BIN="${REPO_DIR}/examples/ttauto"

if [[ ! -x "${TTAUTO_BIN}" ]]; then
  echo "Missing executable: ${TTAUTO_BIN}" >&2
  exit 1
fi

if [[ ! -f "${SCAN_SCRIPT}" || ! -f "${BASELINE_MD}" ]]; then
  echo "Missing scan script or baseline markdown." >&2
  exit 1
fi

scratch="$(mktemp -d)"
trap 'rm -rf "${scratch}"' EXIT

generated_md="${scratch}/ttauto_scan_strata.md"
baseline_before="$(mktemp)"
baseline_after="$(mktemp)"
trap 'rm -rf "${scratch}"; rm -f "${baseline_before}" "${baseline_after}"' EXIT

cp "${BASELINE_MD}" "${baseline_before}"

(
  cd "${scratch}"
  TTAUTO_BIN="${TTAUTO_BIN}" \
  OUTPUT_MD="${generated_md}" \
  SHOW_PROGRESS=0 \
  "${SCAN_SCRIPT}" >/dev/null

  rm -f ./*.m
)

cp "${BASELINE_MD}" "${baseline_after}"

if ! diff -u "${baseline_before}" "${baseline_after}" >/dev/null; then
  echo "Baseline markdown changed during test." >&2
  exit 1
fi

if ! diff -u "${BASELINE_MD}" "${generated_md}" >/dev/null; then
  echo "Generated markdown differs from baseline." >&2
  exit 1
fi
