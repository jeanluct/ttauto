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

# examples/ttauto_strata_braids.md is generated in full by
# examples/ttbraid_strata, prose and all, so regenerating it must reproduce
# it byte for byte.  That guards the braid of every stratum minimiser up to
# seven punctures, and with it the comparison against the published table.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

BASELINE_MD="${REPO_DIR}/examples/ttauto_strata_braids.md"
SWEEP_BIN="${REPO_DIR}/examples/ttbraid_strata"

if [[ ! -x "${SWEEP_BIN}" ]]; then
  echo "Missing executable: ${SWEEP_BIN}" >&2
  exit 1
fi

if [[ ! -f "${BASELINE_MD}" ]]; then
  echo "Missing baseline markdown: ${BASELINE_MD}" >&2
  exit 1
fi

scratch="$(mktemp -d)"
trap 'rm -rf "${scratch}"' EXIT

generated_md="${scratch}/ttauto_strata_braids.md"

(
  cd "${scratch}"
  "${SWEEP_BIN}" > "${generated_md}"
  rm -f ./*.m
)

if ! diff -u "${BASELINE_MD}" "${generated_md}"; then
  echo "Generated markdown differs from baseline." >&2
  exit 1
fi

echo "test_strata_braids_markdown: regenerated output matches the baseline."
